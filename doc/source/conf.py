# -*- coding: utf-8 -*-
#
# PLAMS documentation build configuration file, created by
# sphinx-quickstart2 on Mon Aug 11 16:40:00 2014.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys
from datetime import date

from docutils.parsers.rst.directives.admonitions import BaseAdmonition
from docutils import nodes
from sphinx.util.compat import make_admonition

class Technical(BaseAdmonition):
    node_class = nodes.admonition
    def run(self):
        self.options['class'] = ['technical']
        return make_admonition(
            nodes.admonition, self.name, ["Technical"], self.options,
            self.content, self.lineno, self.content_offset, self.block_text,
            self.state, self.state_machine)

class ADFSuite(BaseAdmonition):
    node_class = nodes.admonition
    def run(self):
        self.options['class'] = ['adfsuite']
        return make_admonition(
            nodes.admonition, self.name, ["ADF Suite"], self.options,
            self.content, self.lineno, self.content_offset, self.block_text,
            self.state, self.state_machine)

def modify_signature(app, what, name, obj, options, signature,
                           return_annotation):
    if signature:
        signature = signature.replace("=u'","='")
    return signature, return_annotation

def setup(app):
    app.add_stylesheet('boxes.css')
    app.add_directive('technical', Technical)
    app.add_directive('adfsuite', ADFSuite)
    app.add_javascript('copybutton.js')
    app.connect('autodoc-process-signature', modify_signature)

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#sys.path.insert(0, os.path.abspath('.'))

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.pngmath',
    'sphinx.ext.viewcode'
]

autodoc_default_flags = ['members', 'private-members', 'special-members']
autodoc_member_order = 'bysource'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'PLAMS'
copyright = u'%i, Scientific Computing & Modelling'%(date.today().year)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '1.0'
# The full version, including alpha/beta/rc tags.
release = 'beta'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []

# The reST default role (used for this markup: `text`) to use for all
# documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
#keep_warnings = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'classic'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "collapsiblesidebar": "true",
    "externalrefs": "true"
}

# Add any paths that contain custom themes here, relative to this directory.
#html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = 'PLAMS documentation'

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = 'PLAMS'

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = '_static/scm_logo_compact.png'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = '_static/favicon.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
#html_extra_path = []

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = False

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'PLAMSdoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
  ('index', 'PLAMS.tex', u'PLAMS Documentation',
   u'Scientific Computing \& Modelling', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'plams', u'PLAMS Documentation',
     [u'Scientific Computing & Modelling'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'PLAMS', u'PLAMS Documentation',
   u'Scientific Computing & Modelling', 'PLAMS', 'One line description of project.',
   'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#texinfo_no_detailmenu = False


# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'python2': ('http://docs.python.org/2.7', None)}

rst_epilog = """
.. |init| replace:: :func:`~scm.plams.common.init`
.. |log| replace:: :func:`~scm.plams.common.log`
.. |load| replace:: :func:`~scm.plams.common.load`
.. |load_all| replace:: :func:`~scm.plams.common.load_all`
.. |finish| replace:: :func:`~scm.plams.common.finish`
.. |add_to_class| replace:: :func:`~scm.plams.common.add_to_class`
.. |add_to_instance| replace:: :func:`~scm.plams.common.add_to_instance`

.. |PlamsError| replace:: :exc:`~scm.plams.errors.PlamsError`
.. |FileError| replace:: :exc:`~scm.plams.errors.FileError`
.. |ResultsError| replace:: :exc:`~scm.plams.errors.ResultsError`
.. |PTError| replace:: :exc:`~scm.plams.errors.PTError`
.. |UnitsError| replace:: :exc:`~scm.plams.errors.UnitsError`
.. |MoleculeError| replace:: :exc:`~scm.plams.errors.MoleculeError`

.. |Job| replace:: :class:`~scm.plams.basejob.Job`
.. |SingleJob| replace:: :class:`~scm.plams.basejob.SingleJob`
.. |MultiJob| replace:: :class:`~scm.plams.basejob.MultiJob`
.. |run| replace:: :meth:`~scm.plams.basejob.Job.run`
.. |prerun| replace:: :meth:`~scm.plams.basejob.Job.prerun`
.. |postrun| replace:: :meth:`~scm.plams.basejob.Job.postrun`

.. |Atom| replace:: :class:`~scm.plams.basemol.Atom`
.. |Bond| replace:: :class:`~scm.plams.basemol.Bond`
.. |Molecule| replace:: :class:`~scm.plams.basemol.Molecule`
.. |BigMolecule| replace:: :class:`~scm.plams.bigmol.BigMolecule`

.. |PeriodicTable| replace:: :class:`~scm.plams.utils.PeriodicTable`
.. |Units| replace:: :class:`~scm.plams.utils.Units`

.. |JobManager| replace:: :class:`~scm.plams.jobmanager.JobManager`
.. |load_job| replace:: :meth:`~scm.plams.jobmanager.JobManager.load_job`

.. |JobRunner| replace:: :class:`~scm.plams.jobrunner.JobRunner`
.. |GridRunner| replace:: :class:`~scm.plams.jobrunner.GridRunner`

.. |Settings| replace:: :class:`~scm.plams.settings.Settings`
.. |Results| replace:: :class:`~scm.plams.results.Results`
.. |KFReader| replace:: :class:`~scm.plams.kftools.KFReader`
.. |KFFile| replace:: :class:`~scm.plams.kftools.KFFile`

.. |SCMJob| replace:: :class:`~scm.plams.scmjob.SCMJob`
.. |SCMResults| replace:: :class:`~scm.plams.scmjob.SCMResults`
.. |ADFJob| replace:: :class:`ADFJob<scm.plams.scmjob.ADFJob>`
.. |ADFResults| replace:: :class:`ADFResults<scm.plams.scmjob.ADFResults>`
.. |BANDJob| replace:: :class:`BANDJob<scm.plams.scmjob.BANDJob>`
.. |BANDResults| replace:: :class:`BANDResults<scm.plams.scmjob.BANDResults>`
.. |DFTBJob| replace:: :class:`DFTBJob<scm.plams.scmjob.DFTBJob>`
.. |DFTBResults| replace:: :class:`DFTBResults<scm.plams.scmjob.DFTBResults>`

.. |DiracJob| replace:: :class:`~scm.plams.diracjob.DiracJob`
.. |DiracResults| replace:: :class:`~scm.plams.diracjob.DiracResults`

.. |RPM| replace:: :ref:`rerun-prevention`
.. |cleaning| replace:: :ref:`cleaning`
.. |pickling| replace:: :ref:`pickling`
.. |binding_decorators| replace:: :ref:`binding-decorators`
.. |parallel| replace:: :ref:`parallel`

"""


