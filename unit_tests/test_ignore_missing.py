from threading import Thread

from scm.plams import Settings, ReentranceError
from scm.plams.core.settings import SupressMissing

glob_var = {
    0: tuple(),
    1: set(),
    2: set(),
    3: set(),
    4: set()
}


class _SupressMissing(SupressMissing):
    @property
    def id(self):
        return id(self.obj.__missing__)

    def __enter__(self):
        """Enter the context manager: modify :meth:`.Settings.__missing__` at the class level."""
        global glob_var
        with self.lock:
            glob_var[1].add(self.id)
            setattr(self.obj, '__missing__', self.missing_new)
            glob_var[2].add(self.id)

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context manager: restore :meth:`.Settings.__missing__` at the class level."""
        global glob_var
        with self.lock:
            glob_var[3].add(self.id)
            setattr(self.obj, '__missing__', self.missing)
            glob_var[4].add(self.id)


class _Settings(Settings):
    ...


_Settings._supress_missing = _SupressMissing(_Settings)
glob_var[0] = (id(_Settings.__missing__),)


def _supress_missing():
    s = _Settings()
    with s.supress_missing():
        pass
    with s.supress_missing():
        pass
    with s.supress_missing():
        pass
    with s.supress_missing():
        pass


def test_supress_missing():
    """Tests for :meth:`Settings.supress_missing`."""
    try:
        thread_list = [Thread(target=_supress_missing) for _ in range(100)]
        for t in thread_list:
            t.start()
        for t in thread_list:
            t.join()

        global glob_var
        for v in glob_var.values():
            assert len(v) == 1

        ref_id = glob_var[0][0]
        assert ref_id in glob_var[1]  # The ID of the unaltered __missing__ method
        assert ref_id not in glob_var[2]  # The ID of the altered __missing__ method
        assert ref_id not in glob_var[3]  # The ID of the altered __missing__ method
        assert ref_id in glob_var[4]  # The ID of the unaltered __missing__ method

    finally:
        for i in range(1, 5):
            glob_var[i] = set()


def test_reentrance():
    """Test that :meth:`Settings.supress_missing` cannot be entered in a reentrant manner."""
    s = Settings()

    try:
        with s.supress_missing():
            with s.supress_missing():
                pass
    except ReentranceError:
        pass
    else:
        raise AssertionError("Failed to raise a 'ReentranceError'")
