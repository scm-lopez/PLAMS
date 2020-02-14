from ..core.settings import Settings
from ..interfaces.adfsuite.ams import AMSResults, AMSJob
from ..interfaces.adfsuite.amsworker import AMSWorkerResults
from ..mol.molecule import Molecule
from ..core.functions import log
from ..core.basejob import MultiJob
import numpy as np

__all__ = ['AMSComparator']


class AMSComparator:
    """A class for running multiple AMS jobs using different settings and comparing the results.

    Initialize an instance with a dictionary of settings, a list of molecules, and a list of properties.

    Call run(), and later access the results with get_results_dict().

    Write to .txt with savetxt()

    Calculate rmsd with get_rmsd()

    add a new column postrun with add_delta_result(), for storing/printing differences in calculated properties

    """

    def __init__(self, settings, molecules, properties="energy"):
        """
            *   ``settings`` -- a dictionary with different AMSJob settings.
            *   ``molecules`` -- a list of molecules. Ideally, each molecule should have a mol.properties.name.
            *   ``properties`` -- a list of properties. For example, properties=['energy', 'gradients'] will call the get_energy() and get_gradients() methods on the AMSResults or AMSWorkerResults after the jobs have been run.

        """
        self.settings, self.molecules = self.ensure_settings_dict_and_molecules_list(settings, molecules)
        #self.molecules = molecules
        self.properties = properties
        if isinstance(self.properties, str):
            self.properties = [self.properties]
        self.all_results_dict = {}
        assert(self.settings is not None)
        assert(self.molecules is not None)
        assert(self.properties is not None)

    def get_results(self):
        """
            If energy, gradients, stresstensor, and other_property were given as properties, then this will
            return a dictionary of the form
            {
                'settings_name_1': {
                    'jobs': [list of AMSJobs],
                    'molecules': [list of PLAMS Molecules],
                    'energy': numpy array of energies,
                    'gradients': numpy array of gradients,
                    'stresstensor': numpy array of stress tensors,
                    'other_property ': ...
                },
                'settings_name_2': {
                    'jobs': [list of AMSJobs],
                    'molecules': [list of PLAMS Molecules],
                    'energy': numpy array of energies,
                    'gradients': numpy array of gradients,
                    'stresstensor': numpy array of stress tensors,
                    'other_property ': ...
                },
            }
            'jobs' is present only if run_mode was 'normal'
        """
        return self.all_results_dict

    def add_delta_result(self, sett1, sett2, new_name=""):
        if not new_name:
            new_name = sett2+"_minus_"+sett1
        assert(new_name not in self.all_results_dict)
        assert(sett1 in self.all_results_dict)
        assert(sett2 in self.all_results_dict)

        self.all_results_dict[new_name]  = dict()
        for prop in self.all_results_dict[sett1]:
            if isinstance(self.all_results_dict[sett1][prop], np.ndarray):
                delta = self.all_results_dict[sett2][prop] - self.all_results_dict[sett1][prop]
                self.all_results_dict[new_name][prop] = delta

    def get_delta(self, sett1, sett2, prop):
        """ Calculates differences between two properties for different settings. 
            sett1: string
            sett2: string
            prop: string
        """
        assert(sett1 in self.all_results_dict)
        assert(sett2 in self.all_results_dict)
        assert(prop in self.all_results_dict[sett1])
        assert(prop in self.all_results_dict[sett2])
        delta = self.all_results_dict[sett2][prop] - self.all_results_dict[sett1][prop]
        return delta

    def get_msd(self, sett1, sett2, prop, per_molecule=True):
        """ Calculate mean squared deviation between two properties for different settings.
        If per_molecule=False, returns a float.
        If per_molecule=True, an array is returned with one entry per molecule """
        assert(sett1 in self.all_results_dict)
        assert(sett2 in self.all_results_dict)
        assert(prop in self.all_results_dict[sett1])
        assert(prop in self.all_results_dict[sett2])
        my_mean = None
        delta = self.get_delta(sett1, sett2, prop)
        squared = delta ** 2

        # if the molecules have different sizes, then the dtype will be object
        if str(squared.dtype) == 'object':
            if per_molecule:
                mean_list = [np.mean(x) for x in squared]
                my_mean = np.array(mean_list)
            else:
                if len(squared[0].shape) > 1: #non-scalar values
                    raveled = np.concatenate([x.ravel() for x in squared])
                    my_mean = np.mean(raveled)
                else:
                    my_mean = np.mean(squared)
        else:
            # this is faster.
            if per_molecule:
                my_mean = np.mean(squared.reshape(squared.shape[0], -1), axis=1)
            else:
                my_mean = np.mean(squared)

        return my_mean

    def get_rmsd(self, sett1, sett2, prop, per_molecule=True):
        assert(sett1 in self.all_results_dict)
        assert(sett2 in self.all_results_dict)
        assert(prop in self.all_results_dict[sett1])
        assert(prop in self.all_results_dict[sett2])
        return np.sqrt(self.get_msd(sett1, sett2, prop, per_molecule))


    def savetxt(self, fname, prop="energy", ravel=False, **kwargs):

        settings_names_list = []
        reshaped_values_list = []
        data_shapes_list = []
        last_dim = 1
        orig_shapes_list = []

        # the assumption here is that all 2D data have the same shape in the second dimension
        # so for example all gradients have second dimension 3, for all molecules
        # the first dimension neeed not the same
        # last_dim is the second dimension of the final data (whether raveled or not) for each setting
        for settings_name, v in self.all_results_dict.items():
            settings_names_list.append(settings_name)
            if prop in v:
                # get the shape of the first item
                orig_shape = v[prop][0].shape
                last_dim = 1
                if len(orig_shape) > 1 and not ravel:
                    last_dim = orig_shape[-1]
                orig_shapes_list.append([x.shape for x in v[prop]])
                data_shapes_list.append([x.reshape((-1, last_dim)).shape for x in v[prop]]) # this should be identical for all settings
                raveled = np.concatenate([x.ravel() for x in v[prop]])
                reshaped = raveled.reshape((-1, last_dim))
                reshaped_values_list.append(reshaped)

        ##### insert data information as the first four columns (prepared as a single string)
        data_information_list = []
        molecule_ids = []
        counter=-1
        for shape, orig_shape, mol in zip(data_shapes_list[0], orig_shapes_list[0], self.molecules):
            counter+=1
            #print(shape)
            if len(orig_shape) > 1 and orig_shape[1] > 1 and last_dim == 1:
                data_index = []
                for outer in range(orig_shape[1]):
                    data_index.extend('{}'.format(x) for x in range(orig_shape[0]))
                data_subindex = ['{}'.format(x) for x in range(orig_shape[1])]*orig_shape[0]
            else:
                data_index = ['{}'.format(x) for x in range(shape[0])]
                data_subindex = ['0']*shape[0]
            thismolecule_name = [mol.properties.name.replace(" ", "_")]*shape[0]


            for i in range(shape[0]):
                complete_string_id = '{} {} {} {}'.format(thismolecule_name[i], counter, data_index[i], data_subindex[i])
                data_information_list.extend([complete_string_id])


        #### create header 
        header = "#mol_name mol_id ind_1(atom) ind_2(component)"
        for s in settings_names_list:
            header += (' '+s)*last_dim

        final_arr = np.concatenate(reshaped_values_list, axis=1)

        fmt = kwargs.get('fmt', '.6e')

        def get_row_string(row_info, row_data):
            ret = str(row_info)+' ' + ' '.join( ('{:'+fmt+'}').format(x) for x in row_data) + '\n'
            return ret

        #def get_columnar_widths():
        #    lines = [header.split()]
        #    for row_info, row_data in zip(data_information_list, final_arr):
        #        lines.extend( [get_row_string(row_info,row_data).split()])
        #    widths = [max(map(len, map(str, col))) for col in zip(*lines)]
        #    return lines, widths

        with open(fname, "w") as f:
            #if columnar:
            #    lines, widths = get_columnar_widths()
            #    for line in lines:
            #        f.write(" ".join((val.ljust(width) for val, width in zip(lines, widths))) + '\n')
            #else:
            f.write(header+'\n')
            for row_info, row_data in zip(data_information_list, final_arr):
                f.write(get_row_string(row_info, row_data))
                    #f.write(str(row_info)+ ' ')
                    ##f.write(' '.join(map(str, row_data))+'\n')
                    #f.write(' '.join( ('{:'+fmt+'}').format(x) for x in row_data))
                    #f.write('\n')


        # np.savetxt works but is very awkawrd with setting the format of the implicit conversion to stirng
        #reshaped_values_list.insert(0, np.array([[molecule_names]]).reshape((-1,1)))
        #np.savetxt(fname, final_arr, header=header, comments='', fmt='%s', **kwargs)
                

    def run(self, run_mode='normal', num_workers=1):
        """
            *   ``run_mode`` -- either 'normal' (leaves files on disk) or 'pipe' (spawns an AMSWorkerPool)
            *   ``num_workers`` -- only if run_mode == 'pipe': the number of workers passed to AMSWorkerPool.
        """

        my_settings_dict, molecules_list = self.settings, self.molecules

        assert(isinstance(self.properties, list))
        assert(isinstance(my_settings_dict, dict) and not isinstance(my_settings_dict, Settings))
        assert(isinstance(molecules_list, list))

        all_results_dict = dict()
        if run_mode == 'normal':
            my_properties = self.clean_properties_list(self.properties, AMSResults) #remove invalid requests
            multijob_list = self.run_multiple_AMS_settings_and_molecules(my_settings_dict, molecules_list)

            results = dict()
            #for engine_name in jobs:
            for multijob in multijob_list:

                engine_name = multijob.name
                results[engine_name] = {}

                for prop in my_properties:
                    results[engine_name][prop] = []

                #for job in jobs[engine_name]:
                for job in multijob.children:
                    for prop in my_properties:
                        p = getattr(job.results, "get_"+prop)()
                        results[engine_name][prop].append(p)

                all_results_dict[engine_name] = dict()
                #all_results_dict[engine_name]['jobs'] = jobs[engine_name]
                #all_results_dict[engine_name]['molecules'] = [ job.results.get_main_molecule() for job in jobs[engine_name] ]
                all_results_dict[engine_name]['jobs'] = multijob.children
                all_results_dict[engine_name]['molecules'] = [ job.results.get_main_molecule() for job in multijob.children ]
                for i, mol in enumerate(molecules_list):
                    all_results_dict[engine_name]['molecules'][i].properties.name = molecules_list[i].properties.name

                for prop in my_properties:
                    all_results_dict[engine_name][prop] = np.array(results[engine_name][prop])

        elif run_mode == 'pipe':
            my_properties = self.clean_properties_list(properties, AMSWorkerResults)
            prop_dict = { x: True for x in my_properties if x != 'energy'}
            results = dict()
            for engine_name, s in my_settings_dict.items():
                if ig('task') in s.input.ams:
                    del s.input.ams[ig('task')]
                all_results_dict[engine_name] = dict()
                all_results_dict[engine_name]['molecules'] = list()
                results[engine_name] = dict()
                for prop in my_properties:
                    results[engine_name][prop] = []
                with AMSWorkerPool(settings=s, num_workers=num_workers) as pool:
                    pool_list = [(mol.properties.name, mol, prop_dict) for mol in molecules_list]
                    pool_results = pool.SinglePoints(pool_list)
                    for result in pool_results:
                        for prop in my_properties:
                            p = getattr(result, "get_"+prop)()
                            results[engine_name][prop].append(p)
                        all_results_dict[engine_name]['molecules'].append(result.get_main_molecule())

                    for i, mol in enumerate(molecules_list):
                        all_results_dict[engine_name]['molecules'][i].properties.name = molecules_list[i].properties.name

                    for prop in my_properties:
                        all_results_dict[engine_name][prop] = np.array(results[engine_name][prop])


        self.all_results_dict = all_results_dict

        return self.all_results_dict


    def ensure_settings_dict_and_molecules_list(self, settings, molecules):
        my_settings_dict = dict()
        if isinstance(settings, Settings):
            if hasattr(settings, 'name'):
                my_settings_dict[s.name] = settings
            else:
                my_settings_dict['unnamed'] = settings
        elif isinstance(settings, list):
            for (i, s) in enumerate(settings):
                if hasattr(s, 'name'):
                    my_settings_dict[str(i)+'_'+s.name] = s
                else:
                    my_settings_dict[str(i)] = s
        elif isinstance(settings, dict) and not isinstance(settings, Settings):
            my_settings_dict = settings.copy()
        else:
            raise RuntimeError("settings must be a Settings, a list of Settings, or a dict['name'] = Settings")
                
        if isinstance(molecules, Molecule):
            molecules = [molecules]
        if isinstance(molecules, list):
            #molecules_dict = dict()
            for i, mol in enumerate(molecules):
                #if hasattr(mol, 'name'):
                    #mol.name = '{}_{}'.format(i, mol.name)
                #else:
                if not hasattr(mol.properties, 'name'):
                    mol.properties.name = '{}_{}'.format(i, mol.get_formula())
                #if hasattr(mol, 'name'):
                #    molecules_dict['{}_{}'.format(i, mol.name)] = mol
                #else:
                #    molecules_dict['{}_{}'.format(i, mol.get_formula())] = mol
            molecules_list = molecules.copy()
        elif isinstance(molecules, dict):
            #molecules_dict = molecules
            for k,v in molecules.items():
                v.name = k
            molecules_list = [v for k,v in molecules.items()]
        else:
            raise RuntimeError("run_comparisons: molecules argument must be one of Molecule, list(Molecule), or dict['mol_name'] => Molecule")

        return my_settings_dict, molecules_list


    def run_multiple_AMS_settings_and_molecules(self, settings, molecules):

        #my_settings_dict, molecules_dict = ensure_settings_dict_and_molecules_dict(settings, molecules)
        #my_settings_dict, molecules_list = self.ensure_settings_dict_and_molecules_list(settings, molecules)
        my_settings_dict, molecules_list = self.settings, self.molecules

        #jobs = dict()
        multijob_list = []
        for engine_name, s in my_settings_dict.items():
            counter = -1
            #jobs[engine_name] = []
            multijob = MultiJob(name=engine_name)
            #for mol_name, mol in molecules_dict.items():
            for mol in molecules_list:
                counter += 1
                #jobname = '{:06d}_{}'.format(counter, mol.properties.name)
                jobname = '{}_{:06d}_{}'.format(engine_name, counter, mol.properties.name) #need the engine name in there to make the job name unique across multijobs
                multijob.children.append( AMSJob(settings=s, molecule=mol, name=jobname) )
                #jobname = '{}_{:06d}_{}'.format(engine_name, counter, mol.properties.name)
                #jobs[engine_name].append( AMSJob(settings=s, molecule=mol, name=jobname) )
                #jobs[engine_name][-1].run()
            multijob_list.append(multijob)
            multijob.run()

        return multijob_list
        #return jobs


    def clean_properties_list(self, properties, cls):
        new_properties = []
        for prop in properties:
            if hasattr(cls, 'get_'+prop):
                new_properties.append(prop)
            else:
                log("Warning, {} does not have a method get_{}. Removing this from the list of requested properties.".format(cls, prop))

        return new_properties
        




