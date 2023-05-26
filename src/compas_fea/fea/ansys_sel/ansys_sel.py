# Author(s): Andrew Liew (github.com/andrewliew), Marius Weber (IBK, ETHZ)

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
from sys import exit

import shutil
import subprocess





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
    
    # Check files exist
    if os.path.exists(check_run_path):
        os.remove(check_run_path)
    else:
        pass

    err_File_ansys=os.path.join(path, name + '.err')        
    if os.path.exists(err_File_ansys):
            os.remove(err_File_ansys)
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

    # Warten bis Ansys Solution finished oder ein Error in Ansys auftritt
    while True:
        
        isFile = os.path.isfile(check_run_path)

        if isFile == True:
            success = True
             #sleep(10)
            break      

        else:            
            test_err=os.path.join(path, name + '.err')                        
            if os.path.isfile(test_err)==True:                  
                with open(test_err) as temp_err:
                    data_err=temp_err.readlines()
                    for line in data_err:
                        if 'ERROR' in line:                                                                                   
                            print('Ansys MAPDL analysis failed - check .err file from Ansys MAPDL')
                            exit()            

    toc = time() - tic
    

    if success:
        if output:
           print('Ansys MAPDL analysis successfull finished in {0:.3f} s'.format(toc))

    else:
        print('Ansys MAPDL analysis failed')
        exit()        

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
        #if steps == 'all':
        #    steps = structure.steps.keys()
        #elif steps == 'last':
        #    steps = [structure.steps_order[-1]]
        #elif type(steps) == str:
        #    steps = [steps]

        for step in steps:
            structure.results[step] = {}
            if structure.steps[step].__name__ == 'GeneralStep':
                
                # Aufbau vektoren fur Node und Element Daten
                elem_infos_list=[]
                rlist = []    
                result_data = []  
                gplist = []


                # extract gerneal element infos
                filename = step + '_elem_infos.txt'

                isfile_filename=str(out_path) + "\\" + filename

                if os.path.isfile(isfile_filename)==True:
                    
                    efile = open(os.path.join(out_path, filename), 'r')    
                    e_i = efile.readlines()                    
                    #leeres Resultat dict.                    
                    elem_infos_dict = {'elem_nr' : {}, 'elem_typ' : {}, 'elem_loc_x_glob_x' : {}, 'elem_loc_x_glob_y' : {},  'elem_loc_x_glob_z' : {}, 'elem_loc_y_glob_x' : {}, 'elem_loc_y_glob_y' : {}, 'elem_loc_y_glob_z' : {}}
                    
                    
                    #print
                    for i in range(len(e_i)):
                        e_i_string = e_i[i].split(',')
                        ele = map(float, e_i_string[1:])
                        key = int(ele[0]) - 1                                                    
                        # Speichert Resultate von fuer Elemente im gesamt Resultatverzeichnis (result_data)
                        elem_infos_dict['elem_nr'][key]=float(ele[0])
                        elem_infos_dict['elem_typ'][key]=float(ele[1])
                        elem_infos_dict['elem_loc_x_glob_x'][key]=float(ele[2])
                        elem_infos_dict['elem_loc_x_glob_y'][key]=float(ele[3])
                        elem_infos_dict['elem_loc_x_glob_z'][key]=float(ele[4])
                        elem_infos_dict['elem_loc_y_glob_x'][key]=float(ele[5])
                        elem_infos_dict['elem_loc_y_glob_y'][key]=float(ele[6])
                        elem_infos_dict['elem_loc_y_glob_z'][key]=float(ele[7])

                    elem_infos_list.append(elem_infos_dict)  

                # Displacements at nodes (node fields)
                # ----------------------------------
                if 'u' in fields or 'all' in fields:
                    filename = step + '_displacements.txt'
                    
                    isfile_filename=str(out_path) + "\\" + filename
                           
                         
                    if os.path.isfile(isfile_filename)==True:          
                        dfile = open(os.path.join(out_path, filename), 'r')
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
                # -------------------------------------------------------
                if 'sf' in fields or 'all' in fields:
                    
                    filename = step + '_shell_forces_moments.txt'

                    isfile_filename=str(out_path) + "\\" + filename

                    if os.path.isfile(isfile_filename)==True:
                        
                        dfile = open(os.path.join(out_path, filename), 'r')    
                        
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
                    
                # Principal Stresses at each GP
                # -------------------------------------------------------
                if 's' in fields or 'all' in fields:
                    
                    # Add stresses elment infos to the structure
                    filename = step + '_stresses_elem_infos.txt'
                    
                    isfile_filename=str(out_path) + "\\" + filename
                                        
                         
                    if os.path.isfile(isfile_filename)==True:          
                        psfile = open(os.path.join(out_path, filename), 'r')
                        ps = psfile.readlines()
                                       
                        
                        stress_dict = {'nr': {},'loc_x_glob_x': {}, 'loc_x_glob_y': {}, 'loc_x_glob_z': {}, 'loc_y_glob_x': {} , 'loc_y_glob_y': {} , 'loc_y_glob_z': {} , 'elem_typ': {} } # sig_c1 and sig_c3  top an einem GP
                        for i in range(len(ps)):
                            psstring = ps[i].split(',')
                            stress = map(float, psstring)
                            key = int(stress[0]) - 1                            
                            stress_dict['nr'][key] = float(stress[0])
                            stress_dict['loc_x_glob_x'][key] = float(stress[1])
                            stress_dict['loc_x_glob_y'][key] = float(stress[2])
                            stress_dict['loc_x_glob_z'][key] = float(stress[3])                            
                            stress_dict['loc_y_glob_x'][key] = float(stress[4])                                                        
                            stress_dict['loc_y_glob_y'][key] = float(stress[5])
                            stress_dict['loc_y_glob_z'][key] = float(stress[6])
                            stress_dict['elem_typ'][key] = float(stress[7])
                        
                        gplist.append(stress_dict)  

                    # Add stresses TOP to the structure
                    filename = step + '_stresses_top.txt'
                    
                    isfile_filename=str(out_path) + "\\" + filename
                                        
                         
                    if os.path.isfile(isfile_filename)==True:          
                        psfile = open(os.path.join(out_path, filename), 'r')
                        ps = psfile.readlines()
                                       
                        
                        stress_dict = {'GP_name_top': {},'elem_nr_top': {}, 'sig_x_top': {}, 'sig_y_top': {}, 'tau_xy_top': {} , 'fcc_eff_top': {} , 'coor_intp_layer_x_top': {} , 'coor_intp_layer_y_top': {}, 'coor_intp_layer_z_top': {}} # sig_c1 and sig_c3  top an einem GP
                        for i in range(len(ps)):
                            psstring = ps[i].split(',')
                            stress = map(float, psstring)
                            key = int(stress[0]) - 1                            
                            stress_dict['GP_name_top'][key] = float(stress[0])
                            stress_dict['elem_nr_top'][key] = float(stress[1])
                            stress_dict['sig_x_top'][key] = float(stress[2])
                            stress_dict['sig_y_top'][key] = float(stress[3])                            
                            stress_dict['tau_xy_top'][key] = float(stress[4])                                                        
                            stress_dict['fcc_eff_top'][key] = float(stress[5])
                            stress_dict['coor_intp_layer_x_top'][key] = float(stress[6])
                            stress_dict['coor_intp_layer_y_top'][key] = float(stress[7])
                            stress_dict['coor_intp_layer_z_top'][key] = float(stress[8])
                        
                        gplist.append(stress_dict)  
                    
                    # Add stresses TOP to the structure      
                    filename = step + '_stresses_bot.txt'
                    
                    isfile_filename=str(out_path) + "\\" + filename
                                        
                         
                    if os.path.isfile(isfile_filename)==True:          
                        psfile = open(os.path.join(out_path, filename), 'r')
                        ps = psfile.readlines()
                                       
                        
                        stress_dict = {'GP_name_bot': {},'elem_nr_bot': {}, 'sig_x_bot': {}, 'sig_y_bot': {}, 'tau_xy_bot': {} , 'fcc_eff_bot': {} , 'coor_intp_layer_x_bot': {} , 'coor_intp_layer_y_bot': {}, 'coor_intp_layer_z_bot': {}} # sig_c1 and sig_c3  top an einem GP
                        for i in range(len(ps)):
                            psstring = ps[i].split(',')
                            stress = map(float, psstring)
                            key = int(stress[0]) - 1
                            stress_dict['GP_name_bot'][key] = float(stress[0])
                            stress_dict['elem_nr_bot'][key] = float(stress[1])
                            stress_dict['sig_x_bot'][key] = float(stress[2])
                            stress_dict['sig_y_bot'][key] = float(stress[3])                            
                            stress_dict['tau_xy_bot'][key] = float(stress[4])                                                        
                            stress_dict['fcc_eff_bot'][key] = float(stress[5])
                            stress_dict['coor_intp_layer_x_bot'][key] = float(stress[6])
                            stress_dict['coor_intp_layer_y_bot'][key] = float(stress[7])
                            stress_dict['coor_intp_layer_z_bot'][key] = float(stress[8])                                                                      
                        
                        gplist.append(stress_dict)
                

                # 
                # 
                            
            #  Speichert nodal and element reuslts in die structure.result dict von Compas FEA. damit die Compas FEA Funktion rhino.plot_data() genutzt werden kann
            

                        
            # Nodal results
            if rlist:
                structure.results[step]['nodal'] = {}
                for rdict in rlist:
                    for key in rdict:
                        structure.results[step]['nodal'][key] = rdict[key] 

            if gplist:
                structure.results[step]['GP'] = {}
                for gpdict in gplist:
                    for key in gpdict:
                        structure.results[step]['GP'][key] = gpdict[key]                         

            # element results results
            if result_data:                
                structure.results[step]['element'] = {}
                structure.results[step]['element'].update(result_data[step]['element'])
            
            # element general infos

            if elem_infos_list:
                structure.results[step]['element_info'] = {}
                for elem_infos_dict in elem_infos_list:
                    for key in elem_infos_dict:
                        structure.results[step]['element_info'][key] = elem_infos_dict[key]  


            toc = time() - tic
            print('Saving Ansys MAPDL results to the structure object successful in {0:.3f} s'.format(toc))

        