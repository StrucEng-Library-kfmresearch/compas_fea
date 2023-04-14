from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas_fea.fea.ansys_sel import Writer
from compas.geometry import length_vector

from subprocess import Popen
from subprocess import PIPE

from time import time
from time import sleep

import json
import os

import shutil
import subprocess


# Author(s): Andrew Liew (github.com/andrewliew), Marius Weber (IBK, ETHZ)


__all__ = [
    'input_generate',
    'extract_data',
    'launch_process',
]


node_fields = ['rf', 'rm', 'u', 'ur', 'cf', 'cm']
element_fields = ['sf', 'sm', 'sk', 'se', 's', 'e', 'pe', 'rbfor', 'ctf']

# -------------------------------------------------------------------------
# Generates the APDL (.inp) file based on the strucutre object
# -------------------------------------------------------------------------
def input_generate(structure, fields, output, lstep, sbstep):
    """ Creates the Ansys .inp file from the Structure object.

    Parameters
    ----------
    structure : obj
        The Structure object to read from.
    fields : list
        Data field requests.
    output : bool
        Print terminal output.

    Returns
    -------
    None

    """

    filename = '{0}{1}.inp'.format(structure.path, structure.name)
    print('')
    print('')
    print('Create Ansys MAPDL input file (.inp)')
    print('--------------------------------------------------------')
    tic = time()
    if isinstance(fields, str):
        fields = [fields]

    if 'u' not in fields:
        fields.append('u')

    with Writer(structure=structure, software='ansys_sel', filename=filename, fields=fields) as writer:
        writer.write_heading() # Writes the Heading in the .inp (APDL) file
        writer.write_nodes() # Writes nodes in the .inp (APDL) file
        writer.write_node_sets() # Writes nodes sets in the .inp (APDL) file
        writer.write_elements(structure) # Writes elements in the .inp (APDL) file
        writer.write_element_sets() # Writes element sets in the .inp (APDL) file
        writer.write_materials() # Writes materials in the .inp (APDL) file
        writer.write_boundary_conditions() # Writes boundary conditions in the .inp (APDL) file        
        writer.write_steps() # Writes steps/solver in the .inp (APDL) file
        writer.write_results(structure,fields, lstep, sbstep) # Writes results in the .inp (APDL) file
    if output:
        toc = time() - tic
        print('Ansys MAPDL input file successfull generated in {0:.3f} s'.format(toc))

# -------------------------------------------------------------------------
# Run ANSYS APDL with the generated APDL (.inp) file
# -------------------------------------------------------------------------
def launch_process(structure, exe, cpus, output):
    """ Runs the analysis through Ansys.

    Parameters
    ----------
    structure : obj
        Structure object.
    exe : str
        ansys exe path to bypass defaults. (not used in the current version)
    cpus : int
        Number of CPU cores to use.
    output : bool
        Print terminal output.

    Returns
    -------
    None

    """
    print('')
    print('')    
    print('Run Ansys MAPDL analysis')
    print('--------------------------------------------------------')
    
    name = structure.name
    path = structure.path
    temp = '{0}{1}/'.format(path, name)
 
    # Analyse
    check_run_path=str(path) + "\\" + 'run_ansys_check.txt'
    
    if os.path.exists(check_run_path):
        os.remove(check_run_path)
    else:
        pass
        
    # Set options
    ansys_path = 'C:\\Program Files\\ANSYS Inc\\v221\\ansys\\bin\\winx64\\ANSYS221.exe'
    lic_str = 'ansys'
    cpus = 1
    inp_path = os.path.join(path, name + '.inp')
    work_dir = os.path.join(path, name + '_output')

    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    out_path = os.path.join(work_dir, name + '.out')

    # Call Ansys
    success = False
     
    launch_string = '\"' + ansys_path + '\" -p ' + lic_str + ' -np ' + str(cpus)
    launch_string += ' -dir \"' + path
    launch_string += '\" -j \"' + name + '\" -s read -l en-us -b -i \"'
    launch_string += inp_path + ' \" -o \"' + out_path + '\"'
    tic = time()
    subprocess = Popen(launch_string, stdout=PIPE, stderr=PIPE, cwd=path, shell=True, env=os.environ)
    #subprocess.wait()

    print('Ansys MAPDL analysis is running ... please wait ... ' )

    while True:
        
        isFile = os.path.isfile(check_run_path)
        
        if isFile == True:
            success = True
             #sleep(10)
            break            
    
    toc = time() - tic
    

    if success:
        if output:
           print('Ansys MAPDL analysis successfull finished in {0:.3f} s'.format(toc))

    else:
        print('Ansys MAPDL analysis failed')

# -------------------------------------------------------------------------
# extract the results from ANSYS APDL and save in the structure 
# -------------------------------------------------------------------------
def extract_data(structure, fields, exe, output, return_data, components):
    
    """ Extract data from the txt files

    Parameters
    ----------
    structure : obj
        Structure object.
    fields : list
        Data field requests.
    exe : str
        Abaqus exe path to bypass defaults. (not used in the current version)
    output : bool
        Print terminal output. (not used in the current version)
    return_data : bool
        Return data back into structure.results.
    components : list
        Specific components to extract from the fields data. (not used in the current version)

    Returns
    -------
    None

    """
    print('')
    print('')    
    print('Extract Ansys MAPDL results to the structure object')
    print('--------------------------------------------------------')

    if return_data: 
        tic = time()   
        steps = structure.steps    
        out_path = os.path.join(structure.path, structure.name + '_output')
        if steps == 'all':
            steps = structure.steps.keys()
        elif steps == 'last':
            steps = [structure.steps_order[-1]]
        elif type(steps) == str:
            steps = [steps]

        for step in steps:
            structure.results[step] = {}
            if structure.steps[step].__name__ == 'GeneralStep':
                rlist = []    
                elist = []  

                # Displacements at nodes (node fields)
                if 'u' in fields or 'all' in fields:
                    filename = step + '_displacements.txt'
                    
                    try:
                        dfile = open(os.path.join(out_path, filename), 'r')
                    except(Exception):
                        return None

                    displacements = dfile.readlines()
                    
                    
                    disp_dict = {'ux': {}, 'uy': {}, 'uz': {}, 'um': {}}
                    for disp in displacements:
                        dstring = disp.split(',')           
                        disp = map(float, dstring[1:])
                        key = int(float(dstring[0])) - 1                        
                        #print(key)
                        disp_dict['ux'][key] = disp[0]
                        disp_dict['uy'][key] = disp[1]
                        disp_dict['uz'][key] = disp[2]
                        # disp_dict['um'][key] = disp[4]
                        disp_dict['um'][key] = length_vector([disp[0], disp[1], disp[2]])
        
                    udict = disp_dict  
                    rlist.append(udict)
                
                # Shell forces/moments at element centroid (element fields)
                if 'sf' in fields or 'all' in fields:
                    filename = step + '_shell_forces_moments.txt'
                    try:
                        dfile = open(os.path.join(out_path, filename), 'r')
                    except(Exception):
                        return None

                    #leeres Resultat dict.
                    result_data = {str(step) : {"element" : {"sf1" : {}, "sf2" : {}, "sf3" : {}, "sf4" : {},  "sf5" : {}, "sm1" : {}, "sm2" : {}, "sm3" : {}, "ele_type" : {},  }}}    
    
                    # displacements = dfile.readlines()
                    shell_forces_moments = dfile.readlines()
                    #print(shell_forces_moments)
                    for f_m in shell_forces_moments:
                        fmstring = f_m.split(',')
                        f_m = map(float, fmstring[2:]) # includes the value sf1, sf2, etc                                        
                        key = int(float(fmstring[1])) - 1
                        
                        # Speichert Resultate von fuer Elemente im gesamt Resultatverzeichnis (result_data)
                        result_data[str(step)]["element"]["sf1"].update({key : {"ip1_sp0" : f_m[0]}})
                        result_data[str(step)]["element"]["sf2"].update({key : {"ip1_sp0" : f_m[1]}})
                        result_data[str(step)]["element"]["sf3"].update({key : {"ip1_sp0" : f_m[2]}})
                        result_data[str(step)]["element"]["sf4"].update({key : {"ip1_sp0" : f_m[3]}})
                        result_data[str(step)]["element"]["sf5"].update({key : {"ip1_sp0" : f_m[4]}})
                        result_data[str(step)]["element"]["sm1"].update({key : {"ip1_sp0" : f_m[5]}})
                        result_data[str(step)]["element"]["sm2"].update({key : {"ip1_sp0" : f_m[6]}})
                        result_data[str(step)]["element"]["sm3"].update({key : {"ip1_sp0" : f_m[7]}})
                        result_data[str(step)]["element"]["ele_type"].update({key : {"ip1_sp0" : f_m[8]}})
                    
            
            #  Speichert nodal and element reuslts in die structure.result dict von Compas FEA. damit die Compas FEA Funktion rhino.plot_data() genutzt werden kann
            # Nodal results
            if rlist:
                structure.results[step]['nodal'] = {}
                for rdict in rlist:
                    for key in rdict:
                        structure.results[step]['nodal'][key] = rdict[key] 

            # element results results

            structure.results[step]['element'] = {}
            structure.results[step]['element'].update(result_data[step]['element'])
            
            
            toc = time() - tic
            print('Saving Ansys MAPDL results to the structure object successful in {0:.3f} s'.format(toc))
            
        