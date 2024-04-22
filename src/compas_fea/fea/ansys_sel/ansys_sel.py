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
def launch_process(structure, exe, cpus, output, ansys_version=None):
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
    ansys_version: string
            Ansys version that shoul be used (e.g. '24' for version 2024 (v241))

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
    if ansys_version==None:
        if os.path.exists('C:\\Program Files\\ANSYS Inc\\v251\\ansys\\bin\\winx64\\ANSYS251.exe'):
            ansys_version='25'
        elif os.path.exists('C:\\Program Files\\ANSYS Inc\\v241\\ansys\\bin\\winx64\\ANSYS241.exe'):
            ansys_version='24'
        elif os.path.exists('C:\\Program Files\\ANSYS Inc\\v231\\ansys\\bin\\winx64\\ANSYS231.exe'):
            ansys_version='23'
        elif os.path.exists('C:\\Program Files\\ANSYS Inc\\v221\\ansys\\bin\\winx64\\ANSYS221.exe'):
            ansys_version='22'
        else: 
            raise Exception("No Ansys Version was found. Please define the ansys version you would like to use.")  
    
    ansys_path = 'C:\\Program Files\\ANSYS Inc\\v{}1\\ansys\\bin\\winx64\\ANSYS{}1.exe'.format(ansys_version,ansys_version)
    print('Ansys Version v'+ansys_version+'1 is used.')
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
    error_found=False # initialising the error found flag
    while True:
        
        isFile = os.path.isfile(check_run_path)

        #test if 'run_ansys_check.txt' is there
        if isFile == True:
            success = True
             #sleep(10)
            break      

        else:  
            #test if '.err' file exists and read it          
            test_err=os.path.join(path, name + '.err')                        
            if os.path.isfile(test_err)==True:                  
                with open(test_err) as temp_err:
                    data_err=temp_err.readlines()
                    for line in data_err:
                        if 'ERROR' in line:                                                                                   
                            print('Ansys MAPDL analysis failed - check .err file from Ansys MAPDL')
                            error_found=True
                            #exit()
                            break   #break out of the for loop - looping through the error file lines
                        else:
                            continue # continue searching lines if no error is found in that line
                if error_found:
                    break # Break out of the while loop because an ERROR was found
                    

    toc = time() - tic
    

    if success:
        if output:
           print('Ansys MAPDL analysis successfull finished in {0:.3f} s'.format(toc))

    else:
        print('Ansys MAPDL analysis failed')
        #exit()  
    return   error_found    

# -------------------------------------------------------------------------
# extract the results from ANSYS APDL and save in the structure 
# -------------------------------------------------------------------------
def extract_data(structure, fields, exe, output, return_data, components, error_found=False):
    
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
    error_found: bool
        Flag that defines weather an error occured during the analysis

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
            structure.results[step] = {} #creates an empty dict for each analysis step
            if structure.steps[step].__name__ == 'GeneralStep':
                
                
                # Aufbau vektoren fur Node und Element Daten
                elem_infos_list=[]
                rlist = []    
                result_data = []  
                gplist = []

                #if error occured write Error into the dict 
                if error_found:
                    structure.results[step]='ERROR'
                    print('The error message was successfully saved to the results dict.')
                    #Attention: Error message is printed into all steps dicts not only int the one where the error occured!

                #if no error occured read out results and write to results dict
                else:
                    # extract gerneal element infos                    
                    #filename = step + '_elem_infos.txt'

                    # ELEM_INFOS.txt
                    # -------------------------------
                    filename = 'elem_infos.txt'

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

                    # ELEM_INFOS_2.txt
                    # -------------------------------
                    filename = 'elem_infos_2.txt'

                    isfile_filename=str(out_path) + "\\" + filename
            
                    if os.path.isfile(isfile_filename)==True:
                        
                        efile = open(os.path.join(out_path, filename), 'r')    
                        e_i = efile.readlines()                    
                        #leeres Resultat dict.                    
                        elem_infos_dict = {'elem_nr' : {}, 'elem_typ' : {}, 'psi_1' : {}, 'psi_2' : {},  'psi_3' : {}, 'psi_4' : {}}
                        
                        
                        #print
                        for i in range(len(e_i)):
                            
                            e_i_string = e_i[i].split(',')
                            ele = map(float, e_i_string[1:])
                            key = int(ele[0]) - 1                                                    
                            # Speichert Resultate von fuer Elemente im gesamt Resultatverzeichnis (result_data)
                            elem_infos_dict['elem_nr'][key]=float(ele[0])
                            elem_infos_dict['elem_typ'][key]=float(ele[1])
                            elem_infos_dict['psi_1'][key]=float(ele[2])
                            elem_infos_dict['psi_2'][key]=float(ele[3])
                            elem_infos_dict['psi_3'][key]=float(ele[4])
                            elem_infos_dict['psi_4'][key]=float(ele[5])

                        elem_infos_list.append(elem_infos_dict)  

                    # ELEM_INFOS_6.txt
                    # -------------------------------
                    filename = 'elem_infos_6.txt'

                    isfile_filename=str(out_path) + "\\" + filename
            
                    if os.path.isfile(isfile_filename)==True:
                        
                        efile = open(os.path.join(out_path, filename), 'r')    
                        e_i = efile.readlines()                    
                        #leeres Resultat dict.                    
                        elem_infos_dict = {'elem_nr' : {}, 'elem_typ' : {}, 'dm1' : {}, 'dm2' : {}, 'dm3' : {}, 'dm4' : {}}
                        
                        
                        #print
                        for i in range(len(e_i)):
                            
                            e_i_string = e_i[i].split(',')
                            ele = map(float, e_i_string[1:])
                            key = int(ele[0]) - 1                                                    
                            # Speichert Resultate von fuer Elemente im gesamt Resultatverzeichnis (result_data)
                            elem_infos_dict['elem_nr'][key]=float(ele[0])
                            elem_infos_dict['elem_typ'][key]=float(ele[1])
                            elem_infos_dict['dm1'][key]=float(ele[2])
                            elem_infos_dict['dm2'][key]=float(ele[3])
                            elem_infos_dict['dm3'][key]=float(ele[4])
                            elem_infos_dict['dm4'][key]=float(ele[5])
    

                        elem_infos_list.append(elem_infos_dict)     

                    # ELEM_INFOS_7.txt
                    # -------------------------------
                    filename = 'elem_infos_7.txt'

                    isfile_filename=str(out_path) + "\\" + filename
            
                    if os.path.isfile(isfile_filename)==True:
                        
                        efile = open(os.path.join(out_path, filename), 'r')    
                        e_i = efile.readlines()                    
                        #leeres Resultat dict.                    
                        elem_infos_dict = {'elem_nr' : {}, 'elem_typ' : {}, 'oo' : {}, 'uu' : {}, 'h_shell' : {}, 'LNr_d06_bot' : {}, 'LNr_d06_top' : {}, 'fcc' : {}}
                        
                        
                        #print
                        for i in range(len(e_i)):
                            
                            e_i_string = e_i[i].split(',')
                            ele = map(float, e_i_string[1:])
                            key = int(ele[0]) - 1                                                    
                            # Speichert Resultate von fuer Elemente im gesamt Resultatverzeichnis (result_data)
                            elem_infos_dict['elem_nr'][key]=float(ele[0])
                            elem_infos_dict['elem_typ'][key]=float(ele[1])
                            elem_infos_dict['oo'][key]=float(ele[2])
                            elem_infos_dict['uu'][key]=float(ele[3])
                            elem_infos_dict['h_shell'][key]=float(ele[4])
                            elem_infos_dict['LNr_d06_bot'][key]=float(ele[5])
                            elem_infos_dict['LNr_d06_top'][key]=float(ele[6])
                            elem_infos_dict['fcc'][key]=float(ele[7])                            
    

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
                                        
                            
                            stress_dict = {'nr': {},'loc_x_glob_x': {}, 'loc_x_glob_y': {}, 'loc_x_glob_z': {}, 'loc_y_glob_x': {} , 'loc_y_glob_y': {} , 'loc_y_glob_z': {} , 'elem_typ': {} , 'fcc': {}, 'k_riss': {}  } # sig_c1 and sig_c3  top an einem GP
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
                                stress_dict['fcc'][key] = float(stress[8])
                                stress_dict['k_riss'][key] = float(stress[9])
                            
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

                        
                        # Add stresses BOT to the structure      
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
  

                    # eps_x, eps_x and eps_xy at each GP
                    # -------------------------------------------------------
                    if 'eps' in fields or 'all' in fields:
                        
                        # Add strains eps_x, eps_y, eps_xy elment infos to the structure
                        filename = step + '_strains_elem_infos.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            psfile = open(os.path.join(out_path, filename), 'r')
                            ps = psfile.readlines()
                                        
                            
                            strain_x_y_xy_dict = {'nr': {},'loc_x_glob_x': {}, 'loc_x_glob_y': {}, 'loc_x_glob_z': {}, 'loc_y_glob_x': {} , 'loc_y_glob_y': {} , 'loc_y_glob_z': {} , 'elem_typ': {} } 
                            for i in range(len(ps)):
                                esp_x_y_y_string = ps[i].split(',')
                                strain_x_y_xy = map(float, esp_x_y_y_string)
                                key = int(strain_x_y_xy[0]) - 1                            
                                strain_x_y_xy_dict['nr'][key] = float(strain_x_y_xy[0])
                                strain_x_y_xy_dict['loc_x_glob_x'][key] = float(strain_x_y_xy[1])
                                strain_x_y_xy_dict['loc_x_glob_y'][key] = float(strain_x_y_xy[2])
                                strain_x_y_xy_dict['loc_x_glob_z'][key] = float(strain_x_y_xy[3])                            
                                strain_x_y_xy_dict['loc_y_glob_x'][key] = float(strain_x_y_xy[4])                                                        
                                strain_x_y_xy_dict['loc_y_glob_y'][key] = float(strain_x_y_xy[5])
                                strain_x_y_xy_dict['loc_y_glob_z'][key] = float(strain_x_y_xy[6])
                                strain_x_y_xy_dict['elem_typ'][key] = float(strain_x_y_xy[7])
                            
                            gplist.append(strain_x_y_xy_dict)  

                        # Add strains eps_x, eps_y, eps_xy elment infos 1 to the structure
                        filename = step + '_strains_elem_infos_1.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            psfile = open(os.path.join(out_path, filename), 'r')
                            ps = psfile.readlines()
                                        
                            
                            strain_x_y_xy_dict = {'nr': {},'ecu': {}, 'k_E': {}, 'k_riss': {}, 'fcc': {}} 
                            for i in range(len(ps)):
                                esp_x_y_y_string = ps[i].split(',')
                                strain_x_y_xy = map(float, esp_x_y_y_string)
                                key = int(strain_x_y_xy[0]) - 1        
                                strain_x_y_xy_dict['nr'][key] = float(strain_x_y_xy[0])                    
                                strain_x_y_xy_dict['ecu'][key] = float(strain_x_y_xy[1])
                                strain_x_y_xy_dict['k_E'][key] = float(strain_x_y_xy[2])
                                strain_x_y_xy_dict['k_riss'][key] = float(strain_x_y_xy[3])
                                strain_x_y_xy_dict['fcc'][key] = float(strain_x_y_xy[4])
                                
                            
                            gplist.append(strain_x_y_xy_dict)                              
                        


                        # Add strains TOP to the structure
                        filename = step + '_eps_x_y_xy_top.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            psfile = open(os.path.join(out_path, filename), 'r')
                            ps = psfile.readlines()
                                        
                            
                            strain_x_y_xy_dict = {'GP_name_top': {},'elem_nr_top': {}, 'eps_x_top': {}, 'eps_y_top': {}, 'eps_xy_top': {}, 'eps_bruch': {} , 'coor_intp_layer_x_top': {} , 'coor_intp_layer_y_top': {}, 'coor_intp_layer_z_top': {}} 
                            for i in range(len(ps)):
                                esp_x_y_y_string = ps[i].split(',')
                                strain_x_y_xy = map(float, esp_x_y_y_string)
                                key = int(strain_x_y_xy[0]) - 1                            
                                strain_x_y_xy_dict['GP_name_top'][key] = float(strain_x_y_xy[0])
                                strain_x_y_xy_dict['elem_nr_top'][key] = float(strain_x_y_xy[1])
                                strain_x_y_xy_dict['eps_x_top'][key] = float(strain_x_y_xy[2])
                                strain_x_y_xy_dict['eps_y_top'][key] = float(strain_x_y_xy[3])                            
                                strain_x_y_xy_dict['eps_xy_top'][key] = float(strain_x_y_xy[4])                                                        
                                strain_x_y_xy_dict['eps_bruch'][key] = float(strain_x_y_xy[5])   
                                strain_x_y_xy_dict['coor_intp_layer_x_top'][key] = float(strain_x_y_xy[6])
                                strain_x_y_xy_dict['coor_intp_layer_y_top'][key] = float(strain_x_y_xy[7])
                                strain_x_y_xy_dict['coor_intp_layer_z_top'][key] = float(strain_x_y_xy[8])
                            
                            gplist.append(strain_x_y_xy_dict)  

                        
                        # Add strains BOT to the structure      
                        filename = step + '_eps_x_y_xy_bot.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            psfile = open(os.path.join(out_path, filename), 'r')
                            ps = psfile.readlines()
                                        
                            
                            strain_x_y_xy_dict = {'GP_name_bot': {},'elem_nr_bot': {}, 'eps_x_bot': {}, 'eps_y_bot': {}, 'eps_xy_bot': {}, 'eps_bruch': {}, 'coor_intp_layer_x_bot': {} , 'coor_intp_layer_y_bot': {}, 'coor_intp_layer_z_bot': {}} 
                            for i in range(len(ps)):
                                esp_x_y_y_string = ps[i].split(',')
                                strain_x_y_xy = map(float, esp_x_y_y_string)
                                key = int(strain_x_y_xy[0]) - 1
                                strain_x_y_xy_dict['GP_name_bot'][key] = float(strain_x_y_xy[0])
                                strain_x_y_xy_dict['elem_nr_bot'][key] = float(strain_x_y_xy[1])
                                strain_x_y_xy_dict['eps_x_bot'][key] = float(strain_x_y_xy[2])
                                strain_x_y_xy_dict['eps_y_bot'][key] = float(strain_x_y_xy[3])                            
                                strain_x_y_xy_dict['eps_xy_bot'][key] = float(strain_x_y_xy[4])                                                                                        
                                strain_x_y_xy_dict['eps_bruch'][key] = float(strain_x_y_xy[5])  
                                strain_x_y_xy_dict['coor_intp_layer_x_bot'][key] = float(strain_x_y_xy[6])
                                strain_x_y_xy_dict['coor_intp_layer_y_bot'][key] = float(strain_x_y_xy[7])
                                strain_x_y_xy_dict['coor_intp_layer_z_bot'][key] = float(strain_x_y_xy[8])    
                                
                            gplist.append(strain_x_y_xy_dict)

                    
                    # Principal Strains at each GP
                    # -------------------------------------------------------
                    if 'eps' in fields or 'all' in fields:
                        
                        # Add strains TOP to the structure
                        filename = step + '_strains_top.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            epsfile = open(os.path.join(out_path, filename), 'r')
                            eps = epsfile.readlines()
                                        
                            
                            strain_dict = {'GP_name_top': {},'elem_nr_top': {}, 'eps_1_top': {}, 'eps_3_top': {}, 'coor_intp_layer_x_top': {} , 'coor_intp_layer_y_top': {}, 'coor_intp_layer_z_top': {}} 
                            for i in range(len(eps)):
                                epsstring = eps[i].split(',')
                                strain = map(float, epsstring)
                                key = int(strain[0]) - 1                            
                                strain_dict['GP_name_top'][key] = float(strain[0])
                                strain_dict['elem_nr_top'][key] = float(strain[1])
                                strain_dict['eps_1_top'][key] = float(strain[2])
                                strain_dict['eps_3_top'][key] = float(strain[3])                            
                                strain_dict['coor_intp_layer_x_top'][key] = float(strain[4])
                                strain_dict['coor_intp_layer_y_top'][key] = float(strain[5])
                                strain_dict['coor_intp_layer_z_top'][key] = float(strain[6])
                            
                            gplist.append(strain_dict)  
                        
                        

                        # Add strain BOT to the structure                          
                        filename = step + '_strains_bot.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            epsfile = open(os.path.join(out_path, filename), 'r')
                            eps = epsfile.readlines()
                                        
                            
                            strain_dict = {'GP_name_bot': {},'elem_nr_bot': {}, 'eps_1_bot': {}, 'eps_3_bot': {}, 'coor_intp_layer_x_bot': {} , 'coor_intp_layer_y_bot': {}, 'coor_intp_layer_z_bot': {}} 
                            for i in range(len(eps)):
                                epsstring = eps[i].split(',')
                                strain = map(float, epsstring)
                                key = int(strain[0]) - 1                            
                                strain_dict['GP_name_bot'][key] = float(strain[0])
                                strain_dict['elem_nr_bot'][key] = float(strain[1])
                                strain_dict['eps_1_bot'][key] = float(strain[2])
                                strain_dict['eps_3_bot'][key] = float(strain[3])                            
                                strain_dict['coor_intp_layer_x_bot'][key] = float(strain[4])
                                strain_dict['coor_intp_layer_y_bot'][key] = float(strain[5])
                                strain_dict['coor_intp_layer_z_bot'][key] = float(strain[6])
                            
                            gplist.append(strain_dict)  
                        
                                        
                    # 
                    # eps_x, eps_x and eps_xy at 06 top and bot 
                    # -------------------------------------------------------
                    if 'eps' in fields or 'all' in fields:
                        
         
                        # Add strains TOP to the structure
                        filename = step + '_eps_x_y_xy_06_bot.txt'
                        
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            psfile = open(os.path.join(out_path, filename), 'r')
                            ps = psfile.readlines()
                                        
                            
                            strain_x_y_xy_dict = {'GP_name_06d_bot': {},'elem_nr_06d_bot': {}, 'eps_x_06d_bot': {}, 'eps_y_06d_bot': {}, 'eps_xy_06d_bot': {}, 'layer_06d_bot': {} , 'coor_intp_layer_x_06d_bot': {} , 'coor_intp_layer_y_06d_bot': {}, 'coor_intp_layer_z_06d_bot': {}} 
                            for i in range(len(ps)):
                                esp_x_y_y_string = ps[i].split(',')
                                strain_x_y_xy = map(float, esp_x_y_y_string)
                                key = int(strain_x_y_xy[0]) - 1                            
                                strain_x_y_xy_dict['GP_name_06d_bot'][key] = float(strain_x_y_xy[0])
                                strain_x_y_xy_dict['elem_nr_06d_bot'][key] = float(strain_x_y_xy[1])
                                strain_x_y_xy_dict['eps_x_06d_bot'][key] = float(strain_x_y_xy[2])
                                strain_x_y_xy_dict['eps_y_06d_bot'][key] = float(strain_x_y_xy[3])                            
                                strain_x_y_xy_dict['eps_xy_06d_bot'][key] = float(strain_x_y_xy[4])                                                        
                                strain_x_y_xy_dict['layer_06d_bot'][key] = float(strain_x_y_xy[5])   
                                strain_x_y_xy_dict['coor_intp_layer_x_06d_bot'][key] = float(strain_x_y_xy[6])
                                strain_x_y_xy_dict['coor_intp_layer_y_06d_bot'][key] = float(strain_x_y_xy[7])
                                strain_x_y_xy_dict['coor_intp_layer_z_06d_bot'][key] = float(strain_x_y_xy[8])
                              
                            gplist.append(strain_x_y_xy_dict)  

                        
          
                        # Add strains TOP to the structure
                        filename = step + '_eps_x_y_xy_06_top.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            psfile = open(os.path.join(out_path, filename), 'r')
                            ps = psfile.readlines()
                                        
                            
                            strain_x_y_xy_dict = {'GP_name_06d_top': {},'elem_nr_06d_top': {}, 'eps_x_06d_top': {}, 'eps_y_06d_top': {}, 'eps_xy_06d_top': {}, 'layer_06d_top': {} , 'coor_intp_layer_x_06d_top': {} , 'coor_intp_layer_y_06d_top': {}, 'coor_intp_layer_z_06d_top': {}} 
                            for i in range(len(ps)):
                                esp_x_y_y_string = ps[i].split(',')
                                strain_x_y_xy = map(float, esp_x_y_y_string)
                                key = int(strain_x_y_xy[0]) - 1                            
                                strain_x_y_xy_dict['GP_name_06d_top'][key] = float(strain_x_y_xy[0])
                                strain_x_y_xy_dict['elem_nr_06d_top'][key] = float(strain_x_y_xy[1])
                                strain_x_y_xy_dict['eps_x_06d_top'][key] = float(strain_x_y_xy[2])
                                strain_x_y_xy_dict['eps_y_06d_top'][key] = float(strain_x_y_xy[3])                            
                                strain_x_y_xy_dict['eps_xy_06d_top'][key] = float(strain_x_y_xy[4])                                                        
                                strain_x_y_xy_dict['layer_06d_top'][key] = float(strain_x_y_xy[5])   
                                strain_x_y_xy_dict['coor_intp_layer_x_06d_top'][key] = float(strain_x_y_xy[6])
                                strain_x_y_xy_dict['coor_intp_layer_y_06d_top'][key] = float(strain_x_y_xy[7])
                                strain_x_y_xy_dict['coor_intp_layer_z_06d_top'][key] = float(strain_x_y_xy[8])
                            
                            gplist.append(strain_x_y_xy_dict)  


                    # Steel stresses reinforcement at each GP
                    # -------------------------------------------------------
                    if 'sig_sr' in fields or 'all' in fields:
                                                
                        # Add stresses elment infos to the structure
                        filename = step + '_sig_sr_elem_infos.txt'
                        
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
                        
                        # Add steel stress reinforcement layer 1 to the structure
                        filename = step + '_sig_sr_1L.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            sig_sr_1_file = open(os.path.join(out_path, filename), 'r')
                            sig_sr_1 = sig_sr_1_file.readlines()
                                        
                            
                            sig_sr_1_dict = {'GP_name_1L': {},'elem_nr_1L': {}, 'sig_sr_1L': {}, 'coor_x_sig_sr_1L': {} , 'coor_y_sig_sr_1L': {}, 'coor_z_sig_sr_1L': {}, 'psi_1L': {}, 'fsy_1L': {}, 'fsu_1L': {}} 
                            for i in range(len(sig_sr_1)):
                                sig_sr_1_string = sig_sr_1[i].split(',')
                                sig_sr_1_stress = map(float, sig_sr_1_string)
                                key = int(sig_sr_1_stress[0]) - 1                            
                                sig_sr_1_dict['GP_name_1L'][key] = float(sig_sr_1_stress[0])
                                sig_sr_1_dict['elem_nr_1L'][key] = float(sig_sr_1_stress[1])
                                sig_sr_1_dict['sig_sr_1L'][key] = float(sig_sr_1_stress[2])                                
                                sig_sr_1_dict['coor_x_sig_sr_1L'][key] = float(sig_sr_1_stress[3])
                                sig_sr_1_dict['coor_y_sig_sr_1L'][key] = float(sig_sr_1_stress[4])
                                sig_sr_1_dict['coor_z_sig_sr_1L'][key] = float(sig_sr_1_stress[5])
                                sig_sr_1_dict['psi_1L'][key] = float(sig_sr_1_stress[6])
                                sig_sr_1_dict['fsy_1L'][key] = float(sig_sr_1_stress[7])
                                sig_sr_1_dict['fsu_1L'][key] = float(sig_sr_1_stress[8])
                            
                            gplist.append(sig_sr_1_dict)  
                        
                        # Add steel stress reinforcement layer 2 to the structure
                        filename = step + '_sig_sr_2L.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            sig_sr_2_file = open(os.path.join(out_path, filename), 'r')
                            sig_sr_2 = sig_sr_2_file.readlines()
                                        
                            
                            sig_sr_2_dict = {'GP_name_2L': {},'elem_nr_2L': {}, 'sig_sr_2L': {}, 'coor_x_sig_sr_2L': {} , 'coor_y_sig_sr_2L': {}, 'coor_z_sig_sr_2L': {}, 'psi_2L': {}, 'fsy_2L': {}, 'fsu_2L': {}}  
                            for i in range(len(sig_sr_2)):
                                sig_sr_2_string = sig_sr_2[i].split(',')
                                sig_sr_2_stress = map(float, sig_sr_2_string)
                                key = int(sig_sr_2_stress[0]) - 1                            
                                sig_sr_2_dict['GP_name_2L'][key] = float(sig_sr_2_stress[0])
                                sig_sr_2_dict['elem_nr_2L'][key] = float(sig_sr_2_stress[1])
                                sig_sr_2_dict['sig_sr_2L'][key] = float(sig_sr_2_stress[2])                                                     
                                sig_sr_2_dict['coor_x_sig_sr_2L'][key] = float(sig_sr_2_stress[3])
                                sig_sr_2_dict['coor_y_sig_sr_2L'][key] = float(sig_sr_2_stress[4])
                                sig_sr_2_dict['coor_z_sig_sr_2L'][key] = float(sig_sr_2_stress[5])
                                sig_sr_2_dict['psi_2L'][key] = float(sig_sr_2_stress[6])
                                sig_sr_2_dict['fsy_2L'][key] = float(sig_sr_2_stress[7])
                                sig_sr_2_dict['fsu_2L'][key] = float(sig_sr_2_stress[8])
                            
                            gplist.append(sig_sr_2_dict)  
                        
                        # Add steel stress reinforcement layer 3 to the structure
                        filename = step + '_sig_sr_3L.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            sig_sr_3_file = open(os.path.join(out_path, filename), 'r')
                            sig_sr_3 = sig_sr_3_file.readlines()
                                        
                            
                            sig_sr_3_dict = {'GP_name_3L': {},'elem_nr_3L': {}, 'sig_sr_3L': {}, 'coor_x_sig_sr_3L': {} , 'coor_y_sig_sr_3L': {}, 'coor_z_sig_sr_3L': {}, 'psi_3L': {}, 'fsy_3L': {}, 'fsu_3L': {}}   
                            for i in range(len(sig_sr_3)):
                                sig_sr_3_string = sig_sr_3[i].split(',')
                                sig_sr_3_stress = map(float, sig_sr_3_string)
                                key = int(sig_sr_3_stress[0]) - 1                            
                                sig_sr_3_dict['GP_name_3L'][key] = float(sig_sr_3_stress[0])
                                sig_sr_3_dict['elem_nr_3L'][key] = float(sig_sr_3_stress[1])
                                sig_sr_3_dict['sig_sr_3L'][key] = float(sig_sr_3_stress[2])                                                            
                                sig_sr_3_dict['coor_x_sig_sr_3L'][key] = float(sig_sr_3_stress[3])
                                sig_sr_3_dict['coor_y_sig_sr_3L'][key] = float(sig_sr_3_stress[4])
                                sig_sr_3_dict['coor_z_sig_sr_3L'][key] = float(sig_sr_3_stress[5])
                                sig_sr_3_dict['psi_3L'][key] = float(sig_sr_3_stress[6])                                
                                sig_sr_3_dict['fsy_3L'][key] = float(sig_sr_3_stress[7])
                                sig_sr_3_dict['fsu_3L'][key] = float(sig_sr_3_stress[8])                                
                            
                            gplist.append(sig_sr_3_dict)                                     

                        
                        # Add steel stress reinforcement layer 4 to the structure
                        filename = step + '_sig_sr_4L.txt'
                        
                        isfile_filename=str(out_path) + "\\" + filename
                                            
                            
                        if os.path.isfile(isfile_filename)==True:          
                            sig_sr_4_file = open(os.path.join(out_path, filename), 'r')
                            sig_sr_4 = sig_sr_4_file.readlines()
                                        
                            
                            sig_sr_4_dict = {'GP_name_4L': {},'elem_nr_4L': {}, 'sig_sr_4L': {}, 'coor_x_sig_sr_4L': {} , 'coor_y_sig_sr_4L': {}, 'coor_z_sig_sr_4L': {}, 'psi_4L': {}, 'fsy_4L': {}, 'fsu_4L': {}}  
                            for i in range(len(sig_sr_4)):
                                sig_sr_4_string = sig_sr_4[i].split(',')
                                sig_sr_4_stress = map(float, sig_sr_4_string)
                                key = int(sig_sr_4_stress[0]) - 1                            
                                sig_sr_4_dict['GP_name_4L'][key] = float(sig_sr_4_stress[0])
                                sig_sr_4_dict['elem_nr_4L'][key] = float(sig_sr_4_stress[1])
                                sig_sr_4_dict['sig_sr_4L'][key] = float(sig_sr_4_stress[2])                                
                                sig_sr_4_dict['coor_x_sig_sr_4L'][key] = float(sig_sr_4_stress[3])
                                sig_sr_4_dict['coor_y_sig_sr_4L'][key] = float(sig_sr_4_stress[4])
                                sig_sr_4_dict['coor_z_sig_sr_4L'][key] = float(sig_sr_4_stress[5])
                                sig_sr_4_dict['psi_4L'][key] = float(sig_sr_4_stress[6])
                                sig_sr_4_dict['fsy_4L'][key] = float(sig_sr_4_stress[7])
                                sig_sr_4_dict['fsu_4L'][key] = float(sig_sr_4_stress[8])   
                            
                            gplist.append(sig_sr_4_dict) 
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
            if not error_found:
                print('Saving Ansys MAPDL results to the structure object successful in {0:.3f} s'.format(toc))

        
