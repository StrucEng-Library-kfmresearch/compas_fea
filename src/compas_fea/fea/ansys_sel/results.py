# Author(s): Andrew Liew (github.com/andrewliew), Marius Weber (IBK, ETHZ)

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import re

import os
import json



__all__ = [
    'Results',
]

dofs = ['x', 'y', 'z', 'xx', 'yy', 'zz']


class Results(object):

    def __init__(self):
        pass

    def write_results(self,structure,fields, lstep, sbstep):
        
        sections = self.structure.sections
        properties = self.structure.element_properties

        # Write part Results
        self.write_section('Results')
        self.blank_line()
        self.write_line('/post1')

        # Delete existing (old) output files
        name = structure.name
        path = structure.path
        filename = name + '_extract.txt'
        out_path = os.path.join(path, name + '_output')
        dir = out_path
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        else:
            for f in os.listdir(dir):
                os.remove(os.path.join(dir, f))
    

        # Check divergency
        self.write_section('Check divergency')
        self.write_line('! ----------------')
        self.blank_line()
        self.write_line('set,last,last')
        self.write_line('*GET,Sol_,ACTIVE,0,SOLU,CNVG')  
        self.blank_line()

        
        # ------------------------------------------------------------------
        # Save general Element infos, ... (MAYBE TO REMOVE)
        # ------------------------------------------------------------------
        elements = self.structure.elements
        num_elements = len(elements)
        properties = self.structure.element_properties
        sections = self.structure.sections
        materials = self.structure.materials
        sets = self.structure.sets
        
        self.write_line('*DIM,elem_infos,ARRAY,{0},33'.format(num_elements))   
        fname = 'elem_infos'
        fname_2 = 'elem_infos_2'
        fname_3 = 'elem_infos_3'
        fname_4 = 'elem_infos_4'
        fname_5 = 'elem_infos_5'
        fname_6 = 'elem_infos_6'
        fname_7 = 'elem_infos_7'
        name_ = 'elem_infos'
            
        # Bestimmung der elementspezifischen Werte
        for key in sorted(properties):
            
            # Import properties, sections, material, type, set, ...
            property = properties[key]
            section = sections[property.section]
            material = materials[property.material]                  
            stype = section.__name__
            elset = property.elset
            element_set = sets[property.elset]
            selection = property.elements if property.elements else sets[elset].selection # Zum elemeset zugehorige Elementnummern
            mtype = material.__name__
            
            # define element type 
            if stype == 'ShellSection':
                stype_number=1 # Schalenelement
            else:
                stype_number=0 # MPC or Other

            # loop all Elements, extract local coor, psi, fsy, fsu, fct, Ec, ecu for cmm usermat 
            for ele_num in sorted(selection):                    
                ele_num_adj=ele_num+1

                # Element infos only for shellelements
                if stype_number == 1: 
                    loc_coor_EV_XA= section.loc_coords_EV_XA
                    loc_coor_EV_YA= section.loc_coords_EV_YA
                    loc_coor_EV_ZA= section.loc_coords_EV_ZA

                    # Bewehrungswinkel psi, fsy, fsu nur bei cmm-usermat; Bei Sandwichmodell gehen Werte in Funktion Sandwichmodel mitein
                    if mtype != 'ElasticIsotropic':  
                        psi1=material.psi1
                        psi2=material.psi2
                        psi3=material.psi3
                        psi4=material.psi4                    
                        
                        fsy1=material.fsy1 
                        fsy2=material.fsy2
                        fsy3=material.fsy3
                        fsy4=material.fsy4

                        fsu1=material.fsu1 
                        fsu2=material.fsu2
                        fsu3=material.fsu3
                        fsu4=material.fsu4    

                        ecu=material.ecu                    
                        k_E=material.k_E
                        k_riss=material.k_riss
                        fcc=material.fcc

                        dm1=material.dm1
                        dm2=material.dm2
                        dm3=material.dm3
                        dm4=material.dm4

                        coo=material.oo
                        cuu=material.uu

                        h_shell=section.h_shell                        
                        nr_layers=section.nr_layers 
                        
                        


                    if element_set.type != 'node':   
                        e_x=loc_coor_EV_XA.get('EV_XA',None)
                        e_y=loc_coor_EV_YA.get('EV_YA',None)
                        e_z=loc_coor_EV_ZA.get('EV_ZA',None)
                        if mtype != 'ElasticIsotropic':  
                            psi_1=psi1.get('psi1',None)
                            psi_2=psi2.get('psi2',None)
                            psi_3=psi3.get('psi3',None)
                            psi_4=psi4.get('psi4',None)

                            fsy_1=fsy1.get('fsy1',None)
                            fsy_2=fsy2.get('fsy2',None)
                            fsy_3=fsy3.get('fsy3',None)
                            fsy_4=fsy4.get('fsy4',None)

                            fsu_1=fsu1.get('fsu1',None)
                            fsu_2=fsu2.get('fsu2',None)
                            fsu_3=fsu3.get('fsu3',None)
                            fsu_4=fsu4.get('fsu4',None) 

                            ecu_C=ecu.get('ecu',None)                   
                            k_E_C=k_E.get('k_E',None)      
                            k_riss_C=k_riss.get('k_riss',None)      
                            fcc_C=fcc.get('fcc',None) 

                            dm_1=dm1.get('dm1',None) 
                            dm_2=dm2.get('dm2',None) 
                            dm_3=dm3.get('dm3',None) 
                            dm_4=dm4.get('dm4',None) 

                            c_oo=coo.get('oo',None) 
                            c_uu=coo.get('oo',None)

                            h_shell_=h_shell.get('h_shell',None)  
                            
                            nn_ = nr_layers.get('nn', None)
                            
                            # values for shear verification
                            delta_h_= h_shell_/nn_

                            d06_bot=((h_shell_-c_oo-dm_4/2)+(h_shell_-c_oo-dm_3/2))/2
                            LNr_d06_bot=round(1/delta_h_*d06_bot)


                            d06_top=((h_shell_-c_uu-dm_1/2)+(h_shell_-c_uu-dm_2/2))/2
                            d06_top_strich=h_shell_-d06_top
                            LNr_d06_top=round(1/delta_h_*d06_top_strich)

                        else:
                            psi_1=0# psi1.get('psi1',None)
                            psi_2=0# psi2.get('psi2',None)
                            psi_3=0# psi3.get('psi3',None)
                            psi_4=0# psi4.get('psi4',None)
                            fsy_1=0
                            fsy_2=0
                            fsy_3=0
                            fsy_4=0
                            fsu_1=0
                            fsu_2=0
                            fsu_3=0
                            fsu_4=0
                            fct_C=0
                            k_E_C=0
                            k_riss_C=0
                            fcc_C=0



                    # Save element infos in vector elem_infos
                    self.write_line('elem_infos({0},1)={1} ! Elementyp e.g. shell=1, mpc and other =0'.format(ele_num_adj,stype_number))   # Elementype     
                    self.write_line('elem_infos({0},2)={1} ! Elementnumber'.format(ele_num_adj,ele_num_adj))   # Elementnumber 
                    self.write_line('elem_infos({0},3)={1}! local x-axis direction global x'.format(ele_num_adj,e_x[0]))     
                    self.write_line('elem_infos({0},4)={1}! local x-axis direction global y'.format(ele_num_adj,e_x[1]))      
                    self.write_line('elem_infos({0},5)={1}! loccal x-axis direction global z'.format(ele_num_adj,e_x[2]))    
                    self.write_line('elem_infos({0},6)={1}! local y-axis direction global x'.format(ele_num_adj,e_y[0]))      
                    self.write_line('elem_infos({0},7)={1}! local y-axis direction global y'.format(ele_num_adj,e_y[1])) 
                    self.write_line('elem_infos({0},8)={1}! local y-axis direction global z'.format(ele_num_adj,e_y[2]))  
                    self.write_line('elem_infos({0},9)={1}! psi1'.format(ele_num_adj,psi_1))  
                    self.write_line('elem_infos({0},10)={1}! psi2'.format(ele_num_adj,psi_2)) 
                    self.write_line('elem_infos({0},11)={1}! psi3'.format(ele_num_adj,psi_3)) 
                    self.write_line('elem_infos({0},12)={1}! psi4'.format(ele_num_adj,psi_4))                                                 
                    self.write_line('elem_infos({0},13)={1}! fsy1'.format(ele_num_adj,fsy_1))           
                    self.write_line('elem_infos({0},14)={1}! fsy2'.format(ele_num_adj,fsy_2))           
                    self.write_line('elem_infos({0},15)={1}! fsy3'.format(ele_num_adj,fsy_3))           
                    self.write_line('elem_infos({0},16)={1}! fsy4'.format(ele_num_adj,fsy_4))           
                    self.write_line('elem_infos({0},17)={1}! fsu1'.format(ele_num_adj,fsu_1))           
                    self.write_line('elem_infos({0},18)={1}! fsu2'.format(ele_num_adj,fsu_2))           
                    self.write_line('elem_infos({0},19)={1}! fsu3'.format(ele_num_adj,fsu_3))           
                    self.write_line('elem_infos({0},20)={1}! fsu4'.format(ele_num_adj,fsu_4))           
                    self.write_line('elem_infos({0},21)={1}! ecu'.format(ele_num_adj,ecu_C))
                    self.write_line('elem_infos({0},22)={1}! k_E'.format(ele_num_adj,k_E_C))
                    self.write_line('elem_infos({0},23)={1}! k_riss'.format(ele_num_adj,k_riss_C))
                    self.write_line('elem_infos({0},24)={1}! dm1'.format(ele_num_adj,dm_1))
                    self.write_line('elem_infos({0},25)={1}! dm2'.format(ele_num_adj,dm_2)) 
                    self.write_line('elem_infos({0},26)={1}! dm3'.format(ele_num_adj,dm_3)) 
                    self.write_line('elem_infos({0},27)={1}! dm4'.format(ele_num_adj,dm_4)) 
                    self.write_line('elem_infos({0},28)={1}! oo'.format(ele_num_adj,c_oo)) 
                    self.write_line('elem_infos({0},29)={1}! uu'.format(ele_num_adj,c_uu)) 
                    self.write_line('elem_infos({0},30)={1}! h_shell'.format(ele_num_adj,h_shell_))                     
                    self.write_line('elem_infos({0},31)={1}! LNr_d06_bot'.format(ele_num_adj,LNr_d06_bot)) 
                    self.write_line('elem_infos({0},32)={1}! LNr_d06_top'.format(ele_num_adj,LNr_d06_top)) 
                    self.write_line('elem_infos({0},33)={1}! fcc'.format(ele_num_adj,fcc_C)) 

                
                else:
                    self.write_line('elem_infos({0},1)={1} ! Elementyp e.g. shell=0, mpc and other =1'.format(ele_num_adj,stype_number))   # Elementype     
                    self.write_line('elem_infos({0},2)={1} ! Elementnumber'.format(ele_num_adj,ele_num_adj))   # Elementnumber 
                    self.write_line('elem_infos({0},3)=0'.format(ele_num_adj))
                    self.write_line('elem_infos({0},4)=0'.format(ele_num_adj))
                    self.write_line('elem_infos({0},5)=0'.format(ele_num_adj)) 
                    self.write_line('elem_infos({0},6)=0'.format(ele_num_adj)) 
                    self.write_line('elem_infos({0},7)=0'.format(ele_num_adj)) 
                    self.write_line('elem_infos({0},8)=0'.format(ele_num_adj))

        # Save elem_infos part 1
        cFile = open(os.path.join(path, filename), 'a')                                       
        self.write_line('allsel')
        self.write_line('*get, nelem, elem,, count ')
        self.write_line('*cfopen,' + out_path + '/' + fname + ',txt ')                
        self.write_line('*do,i,1,nelem ')            
        self.write_line('*CFWRITE, elem_infos, elem_infos(i,2), elem_infos(i,1), elem_infos(i,3), elem_infos(i,4), elem_infos(i,5), elem_infos(i,6), elem_infos(i,7), elem_infos(i,8)')
        self.write_line('*Enddo ')
        self.write_line('ESEL, ALL ')
        self.write_line('ETABLE, ERAS ')
        self.write_line('! ')
        cFile.close()

        # Save elem_infos part 2
        cFile = open(os.path.join(path, filename), 'a')
        self.write_line('allsel')
        self.write_line('*get, nelem, elem,, count ')
        self.write_line('*cfopen,' + out_path + '/' + fname_2 + ',txt ')                
        self.write_line('*do,i,1,nelem ')            
        self.write_line('*CFWRITE, elem_infos, elem_infos(i,2), elem_infos(i,1), elem_infos(i,9), elem_infos(i,10), elem_infos(i,11), elem_infos(i,12)')
        self.write_line('*Enddo ')
        self.write_line('ESEL, ALL ')
        self.write_line('ETABLE, ERAS ')
        self.write_line('! ')
        cFile.close() 

        # Save elem_infos part 3
        cFile = open(os.path.join(path, filename), 'a')
        self.write_line('allsel')
        self.write_line('*get, nelem, elem,, count ')
        self.write_line('*cfopen,' + out_path + '/' + fname_3 + ',txt ')                
        self.write_line('*do,i,1,nelem ')            
        self.write_line('*CFWRITE, elem_infos, elem_infos(i,2), elem_infos(i,1), elem_infos(i,13), elem_infos(i,14), elem_infos(i,15), elem_infos(i,16)')
        self.write_line('*Enddo ')
        self.write_line('ESEL, ALL ')
        self.write_line('ETABLE, ERAS ')
        self.write_line('! ')
        cFile.close()    

        # Save elem_infos part 4
        cFile = open(os.path.join(path, filename), 'a')
        self.write_line('allsel')
        self.write_line('*get, nelem, elem,, count ')
        self.write_line('*cfopen,' + out_path + '/' + fname_4 + ',txt ')                
        self.write_line('*do,i,1,nelem ')            
        self.write_line('*CFWRITE, elem_infos, elem_infos(i,2), elem_infos(i,1), elem_infos(i,17), elem_infos(i,18), elem_infos(i,19), elem_infos(i,20)')
        self.write_line('*Enddo ')
        self.write_line('ESEL, ALL ')
        self.write_line('ETABLE, ERAS ')
        self.write_line('! ')
        cFile.close()       


        # Save elem_infos part 5
        cFile = open(os.path.join(path, filename), 'a')
        self.write_line('allsel')
        self.write_line('*get, nelem, elem,, count ')
        self.write_line('*cfopen,' + out_path + '/' + fname_5 + ',txt ')                
        self.write_line('*do,i,1,nelem ')            
        self.write_line('*CFWRITE, elem_infos, elem_infos(i,2), elem_infos(i,21), elem_infos(i,22), elem_infos(i,23), elem_infos(i,24)')
        self.write_line('*Enddo ')
        self.write_line('ESEL, ALL ')
        self.write_line('ETABLE, ERAS ')
        self.write_line('! ')
        cFile.close()

        # Save elem_infos part 6
        cFile = open(os.path.join(path, filename), 'a')
        self.write_line('allsel')
        self.write_line('*get, nelem, elem,, count ')
        self.write_line('*cfopen,' + out_path + '/' + fname_6 + ',txt ')                
        self.write_line('*do,i,1,nelem ')            
        self.write_line('*CFWRITE, elem_infos, elem_infos(i,2), elem_infos(i,1), elem_infos(i,24), elem_infos(i,25), elem_infos(i,26), elem_infos(i,27)')
        self.write_line('*Enddo ')
        self.write_line('ESEL, ALL ')
        self.write_line('ETABLE, ERAS ')
        self.write_line('! ')
        cFile.close()         

        # Save elem_infos part 7
        cFile = open(os.path.join(path, filename), 'a')
        self.write_line('allsel')
        self.write_line('*get, nelem, elem,, count ')
        self.write_line('*cfopen,' + out_path + '/' + fname_7 + ',txt ')                
        self.write_line('*do,i,1,nelem ')            
        self.write_line('*CFWRITE, elem_infos, elem_infos(i,2), elem_infos(i,1), elem_infos(i,28), elem_infos(i,29), elem_infos(i,30), elem_infos(i,31), elem_infos(i,32), elem_infos(i,33)')
        self.write_line('*Enddo ')
        self.write_line('ESEL, ALL ')
        self.write_line('ETABLE, ERAS ')
        self.write_line('! ')
        cFile.close()                                   
                                          
        
        # ------------------------------------------------------------------
        # Schleife uber alle angegebenen lstep (e.g. 'step_3')
        # ------------------------------------------------------------------
        self.write_section('Start extracting selected results for each selected step')
        
        # Define times and steps
        sbstep=sbstep # equal to increments        
        step_list=structure.steps_order        

        for single_lstep in lstep:            
            lstep_index=step_list.index(single_lstep)         
            self.write_line('set,{0},{1}'.format(lstep_index,sbstep))
            print('---------------------')
            
            step_name = structure.steps_order[lstep_index]

            # ------------------------------------------------------------------
            # WRITE DATA AT NODES
            # ------------------------------------------------------------------

            # Write displacements at nodes
            # ------------------------------------------------------------------
            if 'u' in fields or 'all' in fields:
                fname = str(step_name) + '_' + 'displacements'
                name_ = 'nds_d'
                name_x = 'dispX'
                name_y = 'dispY'
                name_z = 'dispZ'

                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write Displacements')
                self.blank_line()
                self.write_line('nsel, all')
                self.write_line('*get,numNodes,node,,count')
                self.write_line('*set,' + name_x + ',')
                self.write_line('*dim,' + name_x + ',array,numNodes,1')
                self.write_line('*set,' + name_y + ',')
                self.write_line('*dim,' + name_y + ',array,numNodes,1')
                self.write_line('*set,' + name_z + ',')
                self.write_line('*dim,' + name_z + ',array,numNodes,1')
                self.write_line('*dim,' + name_ + ', ,numNodes')
                self.write_line('*VGET, ' + name_x + ', node, all, u, X,,,2') # Node displacement in x-direction
                self.write_line('*VGET, ' + name_y + ', node, all, u, Y,,,2') # Node displacement in y-direction
                self.write_line('*VGET, ' + name_z + ', node, all, u, Z,,,2') # Node displacement in z-direction
                self.write_line('!')
                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_x + '(1) , \',\' , ' + name_y + '(1) , \',\' ,' + name_z + '(1)')
                self.write_line('(          F100000.0,       A,       ES,           A,          ES,          A,      ES)')
                self.write_line('*cfclose \n')
                self.write_line('!')
                self.write_line('!')
                cFile.close()

            # ------------------------------------------------------------------
            # WRITE DATA AT SHELL MIDPOINT
            # ------------------------------------------------------------------                

            # Write shell forces and montents at shell mitdpoint (only for Shell elements)
            # ------------------------------------------------------------------
            if 'sf' in fields or 'all' in fields:

                fname = str(step_name) + '_' + 'shell_forces_moments'
                name_ = 'ele_sf'

                # Abfuellen der Resultates
                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write Shell forces and moments')
                self.blank_line()
                self.write_line('allsel')
                self.write_line('ETABLE, ERAS ')
                self.write_line('ETABLE, , SMISC, 1 ')  # N11 In-plane forces (per unit length)
                self.write_line('ETABLE, , SMISC, 2 ')  # N22 In-plane forces (per unit length)
                self.write_line('ETABLE, , SMISC, 3 ')  # N12 In-plane forces (per unit length)

                self.write_line('ETABLE, , SMISC, 4 ')  # M11 Out-of-plane moments (per unit length)
                self.write_line('ETABLE, , SMISC, 5 ')  # M22 Out-of-plane moments (per unit length)
                self.write_line('ETABLE, , SMISC, 6 ')  # M12 Out-of-plane moments (per unit length)

                self.write_line('ETABLE, , SMISC, 7 ')  # Q13 Out-of-plane moments (per unit length)
                self.write_line('ETABLE, , SMISC, 8 ')  # Q23 Out-of-plane moments (per unit length)

                self.write_line('*get, nelem, elem,, count ')
                self.write_line('*dim, enum, array, nelem, 1 ')
                self.write_line('*dim, eforces, array, nelem, 9  ')
                self.write_line('*vget, enum, ELEM, , ELIST ')
            

                self.write_line('*do,i,1,nelem ')
                self.write_line('*if,elem_infos(i,1),EQ,1,THEN  ')             # Abfullen bei Schalenelemente
                self.write_line('*get, eforces(i, 1), ETAB, 1, ELEM, enum(i) ') # N11 In-plane forces (per unit length)
                self.write_line('*get, eforces(i, 2), ETAB, 2, ELEM, enum(i) ') # N22 In-plane forces (per unit length)
                self.write_line('*get, eforces(i, 3), ETAB, 3, ELEM, enum(i) ') # N12 In-plane forces (per unit length)
                self.write_line('*get, eforces(i, 4), ETAB, 7, ELEM, enum(i) ') # Q13 Out-of-plane moments (per unit length)
                self.write_line('*get, eforces(i, 5), ETAB, 8, ELEM, enum(i) ') # Q23 Out-of-plane moments (per unit length)
                self.write_line('*get, eforces(i, 6), ETAB, 4, ELEM, enum(i) ') # M11 Out-of-plane moments (per unit length)
                self.write_line('*get, eforces(i, 7), ETAB, 5, ELEM, enum(i) ') # M22 Out-of-plane moments (per unit length)
                self.write_line('*get, eforces(i, 8), ETAB, 6, ELEM, enum(i) ') # M12 Out-of-plane moments (per unit length)       
                self.write_line('eforces(i, 9)=elem_infos(i,1)') 
                self.write_line('*else')    
                self.write_line('*endif')    
                self.write_line('*Enddo ')

                fname = str(step_name) + '_' + 'shell_forces_moments'
                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt ')
                #self.write_line('*CFWRITE, stresses, E number, Sm11, Sm22, Sm12, Sb11 \n')
                self.write_line('*do,i,1,nelem ')
                #self.write_line('*CFWRITE, stresses, enum(i,1), eforces(i,1), eforces(i,2),eforces(i,3), eforces(i,2) \n')
                self.write_line('*CFWRITE, shell_forces_moments, enum(i,1), eforces(i,1), eforces(i,2),eforces(i,3), eforces(i,4), eforces(i,5), eforces(i,6), eforces(i,7), eforces(i,8), eforces(i,9)   ')
                self.write_line('*Enddo ')
                self.write_line('ESEL, ALL ')
                self.write_line('ETABLE, ERAS ')
                self.write_line('! ')
                cFile.close()
            
           
            # ------------------------------------------------------------------
            # WRITE DATA AT GP
            # ------------------------------------------------------------------    
        
            # Write Stresses at top and bottom Elementlayer at every GP 
            # ------------------------------------------------------------------               
            if 's' in fields or 'all' in fields:

                # Benotigte Element infos abspeichern fur stresses (MAYBE TO REMOVE)
                # --------------------------------------------
                fname = str(step_name) + '_' + 'stresses_elem_infos'
                name_ = 'e_info_GP'
                name_loc_x_glob_x ='loc_x_glob_x'
                name_loc_x_glob_y ='loc_x_glob_y'
                name_loc_x_glob_z ='loc_x_glob_z'
                name_loc_y_glob_x ='loc_y_glob_x'
                name_loc_y_glob_y ='loc_y_glob_y'
                name_loc_y_glob_z ='loc_y_glob_z'
                name_elem_typ ='elem_typ'
                name_fcc ='fcc'
                name_kriss ='k_riss'

                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write element infos for stresses in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,loc_x_glob_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_x_glob_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_x_glob_z,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_y_glob_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_y_glob_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_y_glob_z,ARRAY,4*NrE,1')
                self.write_line('*DIM,elem_typ,ARRAY,4*NrE,1')
                self.write_line('*DIM,fcc,ARRAY,4*NrE,1')
                self.write_line('*DIM,k_riss,ARRAY,4*NrE,1')

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')                             
                            
                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')


                self.write_line('loc_x_glob_x(aux,1)=elem_infos(ii,3)') 
                self.write_line('loc_x_glob_y(aux,1)=elem_infos(ii,4)') 
                self.write_line('loc_x_glob_z(aux,1)=elem_infos(ii,5)') 
                self.write_line('loc_y_glob_x(aux,1)=elem_infos(ii,6)') 
                self.write_line('loc_y_glob_y(aux,1)=elem_infos(ii,7)') 
                self.write_line('loc_y_glob_z(aux,1)=elem_infos(ii,8)') 
                self.write_line('elem_typ(aux,1)=elem_infos(ii,1)') 
                self.write_line('fcc(aux,1)=elem_infos(ii,24)')    
                self.write_line('k_riss(aux,1)=elem_infos(ii,23)')         

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR')
                self.write_line('*else')    
                self.write_line('*endif')  
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')

                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_loc_x_glob_x + '(1) , \',\' , ' + name_loc_x_glob_y + '(1) , \',\' ,' + name_loc_x_glob_z + '(1) , \',\' ,'+ name_loc_y_glob_x + '(1) , \',\' ,' + name_loc_y_glob_y + '(1) , \',\' ,' + name_loc_y_glob_z + '(1) , \',\' ,' + name_elem_typ + '(1) , \',\' ,' + name_fcc + '(1) , \',\' ,'  + name_kriss + '(1) ')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('!')
                self.write_line('!')
                cFile.close()                   


                # Top Layer (nn)
                # --------------------------------------------
                fname = str(step_name) + '_' + 'stresses_top'
                name_ = 's_GP'
                name_elem_nr = 'elem_nr'
                name_sig_c1 = 'sig_c1'
                name_sig_c3 = 'sig_c3'
                name_sig_x = 'sig_x'
                name_sig_y = 'sig_y'
                name_tau_xy = 'tau_xy'
                name_fcc_eff = 'fcc_eff'
                name_coor_intp_layer_x ='coor_intp_layer_x'
                name_coor_intp_layer_y ='coor_intp_layer_y'
                name_coor_intp_layer_z ='coor_intp_layer_z'

        
                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write stresses in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,sig_c1,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_c3,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,tau_xy,ARRAY,4*NrE,1')
                self.write_line('*DIM,fcc_eff,ARRAY,4*NrE,1')

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                self.write_line('*DIM,coor_intp_layer_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_z,ARRAY,4*NrE,1')               
  
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')          
                
                # (NrT: Shell ID, NrL=Numbers of Layer) 
                self.write_line('*GET, NrT, ELEM,N_E(ii), attr, secn') # gibt zu einem Element zugehorgie secnum
                self.write_line('*GET, NrL, SHEL, NrT, Prop,NLAY ') # NrT is equal to secnum (not Element number)
                            

                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')
                self.write_line('*GET,s_GP(aux,1),NODE,N_N(kk),SVAR,1 ')
                

                self.write_line('LAYER,NrL')
                self.write_line('elem_nr(aux,1)=N_E(ii)')
                self.write_line('*GET,usedmodel_check,NODE,N_N(kk),SVAR,1 ')
                
                self.write_line('*IF,usedmodel_check,EQ,1,THEN ')                                
                self.write_line('*GET,sig_c1(aux,1),NODE,N_N(kk),SVAR,8')                
                self.write_line('*GET,sig_c3(aux,1),NODE,N_N(kk),SVAR,9')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68') 
                self.write_line('*GET,fcc_eff(aux,1),NODE,N_N(kk),SVAR,12') 

                self.write_line('*ELSEIF,usedmodel_check,EQ,2,THEN')
                self.write_line('*GET,sig_c1(aux,1),NODE,N_N(kk),SVAR,29 ')
                self.write_line('*GET,sig_c3(aux,1),NODE,N_N(kk),SVAR,30 ')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68')                 
                self.write_line('*GET,fcc_eff(aux,1),NODE,N_N(kk),SVAR,33 ') 

                self.write_line('*ELSEIF,usedmodel_check,EQ,3,THEN') # linear elastisch
                self.write_line('*GET,sig_x,NODE,N_N(kk),SVAR,66')
                self.write_line('*GET,sig_y,NODE,N_N(kk),SVAR,67')
                self.write_line('*GET,sig_xy,NODE,N_N(kk),SVAR,68')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68') 
                self.write_line('sig_c1(aux,1)=(sig_x+sig_y)/2+SQRT((1/2)*((sig_x+sig_y)*(sig_x+sig_y)+(2*sig_xy)*(2*sig_xy)))')
                self.write_line('sig_c3(aux,1)=(sig_x+sig_y)/2-SQRT((1/2)*((sig_x+sig_y)*(sig_x+sig_y)+(2*sig_xy)*(2*sig_xy)))')
                self.write_line('fcc_eff(aux,1)=elem_infos(ii,24)')

                self.write_line('*ELSEIF,usedmodel_check,EQ,4,THEN')
                self.write_line('sig_c1(aux,1)=0')
                self.write_line('*GET,sig_c3(aux,1),NODE,N_N(kk),SVAR,18 ')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68')                
                self.write_line('*GET,fcc_eff(aux,1),NODE,N_N(kk),SVAR,19') 

                self.write_line('*ELSEIF,usedmodel_check,EQ,5,THEN') # Zug Zug
                self.write_line('*GET,sig_c1(aux,1),NODE,N_N(kk),SVAR,39')
                self.write_line('sig_c3(aux,1)=0')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68')                
                self.write_line('fcc_eff(aux,1)=elem_infos(ii,24)')
                
                self.write_line('*ENDIF')

                self.write_line('*GET,coor_intp_layer_x(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_intp_layer_y(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_intp_layer_z(aux,1),NODE,N_N(kk),SVAR,65')               

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR')
                self.write_line('*else')    
                self.write_line('*endif')  
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_x + '(1) , \',\' ,' + name_sig_y + '(1) , \',\' ,'+ name_tau_xy + '(1) , \',\' ,' + name_fcc_eff + '(1) , \',\' ,'+ name_coor_intp_layer_x + '(1) , \',\' ,'  + name_coor_intp_layer_y + '(1) , \',\' ,' + name_coor_intp_layer_z + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('!')
                self.write_line('!')
                cFile.close()      

                # Bot Layer (1)
                # --------------------------------------------
                fname = str(step_name) + '_' + 'stresses_bot'
                name_ = 's_GP'
                name_elem_nr = 'elem_nr'
                name_sig_c1 = 'sig_c1'
                name_sig_c3 = 'sig_c3'
                name_sig_x = 'sig_x'
                name_sig_y = 'sig_y'
                name_tau_xy = 'tau_xy'
                name_fcc_eff = 'fcc_eff'
                name_coor_intp_layer_x ='coor_intp_layer_x'
                name_coor_intp_layer_y ='coor_intp_layer_y'
                name_coor_intp_layer_z ='coor_intp_layer_z'

                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write stresses')
                self.blank_line()


                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,sig_c1,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_c3,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,tau_xy,ARRAY,4*NrE,1')
                self.write_line('*DIM,fcc_eff,ARRAY,4*NrE,1')

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                self.write_line('*DIM,coor_intp_layer_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_z,ARRAY,4*NrE,1')
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE')
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')

                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')                         
                  
                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')
                self.write_line('*GET,s_GP(aux,1),NODE,N_N(kk),SVAR,1 ')
                
                self.write_line('LAYER,1')
                self.write_line('elem_nr(aux,1)=N_E(ii) ')
                self.write_line('*GET,usedmodel_check,NODE,N_N(kk),SVAR,1 ')
                
                self.write_line('*IF,usedmodel_check,EQ,1,THEN ')                
                self.write_line('*GET,sig_c1(aux,1),NODE,N_N(kk),SVAR,8')                
                self.write_line('*GET,sig_c3(aux,1),NODE,N_N(kk),SVAR,9')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68') 
                self.write_line('*GET,fcc_eff(aux,1),NODE,N_N(kk),SVAR,12') 

                self.write_line('*ELSEIF,usedmodel_check,EQ,2,THEN')
                self.write_line('*GET,sig_c1(aux,1),NODE,N_N(kk),SVAR,29 ')
                self.write_line('*GET,sig_c3(aux,1),NODE,N_N(kk),SVAR,30 ')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68')                 
                self.write_line('*GET,fcc_eff(aux,1),NODE,N_N(kk),SVAR,33 ') 

                self.write_line('*ELSEIF,usedmodel_check,EQ,3,THEN') # linear elastisch
                self.write_line('*GET,sig_x,NODE,N_N(kk),SVAR,66')
                self.write_line('*GET,sig_y,NODE,N_N(kk),SVAR,67')
                self.write_line('*GET,sig_xy,NODE,N_N(kk),SVAR,68')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68') 
                self.write_line('sig_c1(aux,1)=(sig_x+sig_y)/2+SQRT((1/2)*((sig_x+sig_y)*(sig_x+sig_y)+(2*sig_xy)*(2*sig_xy)))')
                self.write_line('sig_c3(aux,1)=(sig_x+sig_y)/2-SQRT((1/2)*((sig_x+sig_y)*(sig_x+sig_y)+(2*sig_xy)*(2*sig_xy)))')
                self.write_line('fcc_eff(aux,1)=elem_infos(ii,24)')

                self.write_line('*ELSEIF,usedmodel_check,EQ,4,THEN')
                self.write_line('sig_c1(aux,1)=0')
                self.write_line('*GET,sig_c3(aux,1),NODE,N_N(kk),SVAR,18 ')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68')                
                self.write_line('*GET,fcc_eff(aux,1),NODE,N_N(kk),SVAR,19') 

                self.write_line('*ELSEIF,usedmodel_check,EQ,5,THEN') # Zug Zug
                self.write_line('*GET,sig_c1(aux,1),NODE,N_N(kk),SVAR,39')
                self.write_line('sig_c3(aux,1)=0')
                self.write_line('*GET,sig_x(aux,1),NODE,N_N(kk),SVAR,66')                
                self.write_line('*GET,sig_y(aux,1),NODE,N_N(kk),SVAR,67')                
                self.write_line('*GET,tau_xy(aux,1),NODE,N_N(kk),SVAR,68')                
                self.write_line('fcc_eff(aux,1)=elem_infos(ii,24)')
                
                self.write_line('*ENDIF')

                self.write_line('*GET,coor_intp_layer_x(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_intp_layer_y(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_intp_layer_z(aux,1),NODE,N_N(kk),SVAR,65')

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR')
                self.write_line('*else')    
                self.write_line('*endif')                  
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')

                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_x + '(1) , \',\' ,' + name_sig_y + '(1) , \',\' ,'+ name_tau_xy + '(1) , \',\' ,' + name_fcc_eff + '(1) , \',\' ,'+ name_coor_intp_layer_x + '(1) , \',\' ,'  + name_coor_intp_layer_y + '(1) , \',\' ,' + name_coor_intp_layer_z + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('!')
                self.write_line('!')
                cFile.close()      
          
            else:
                pass 
                              
                


            # Write eps_x, eps_y and eps_xy strains at top and bottom Elementlayer at every GP
            # -------------------------------------------------------------------------------   
            if 'eps' in fields or 'all' in fields:
    
                # Benotigte Element infos abspeichern fur stresses (MAYBE TO REMOVE)
                # --------------------------------------------
                fname = str(step_name) + '_' + 'strains_elem_infos'
                fname_1 = str(step_name) + '_' + 'strains_elem_infos_1'
                name_ = 'eps_info_GP'                
                name_loc_x_glob_x ='loc_x_glob_x'
                name_loc_x_glob_y ='loc_x_glob_y'
                name_loc_x_glob_z ='loc_x_glob_z'
                name_loc_y_glob_x ='loc_y_glob_x'
                name_loc_y_glob_y ='loc_y_glob_y'
                name_loc_y_glob_z ='loc_y_glob_z'
                name_elem_typ ='elem_typ'
                name_ecu ='ecu'
                name_kE ='k_E'
                name_kriss ='k_riss'
                name_fcc ='fcc'
                
                

                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write element infos for eps_x, eps_y and eps_xy in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,loc_x_glob_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_x_glob_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_x_glob_z,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_y_glob_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_y_glob_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_y_glob_z,ARRAY,4*NrE,1')
                self.write_line('*DIM,elem_typ,ARRAY,4*NrE,1')
                self.write_line('*DIM,ecu,ARRAY,4*NrE,1')
                self.write_line('*DIM,k_E,ARRAY,4*NrE,1')
                self.write_line('*DIM,k_riss,ARRAY,4*NrE,1')
                self.write_line('*DIM,fcc,ARRAY,4*NrE,1')
                

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')                             
                            
                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')


                self.write_line('loc_x_glob_x(aux,1)=elem_infos(ii,3)') 
                self.write_line('loc_x_glob_y(aux,1)=elem_infos(ii,4)') 
                self.write_line('loc_x_glob_z(aux,1)=elem_infos(ii,5)') 
                self.write_line('loc_y_glob_x(aux,1)=elem_infos(ii,6)') 
                self.write_line('loc_y_glob_y(aux,1)=elem_infos(ii,7)') 
                self.write_line('loc_y_glob_z(aux,1)=elem_infos(ii,8)') 
                self.write_line('elem_typ(aux,1)=elem_infos(ii,1)') 
                self.write_line('ecu(aux,1)=elem_infos(ii,21)') 
                self.write_line('k_E(aux,1)=elem_infos(ii,22)') 
                self.write_line('k_riss(aux,1)=elem_infos(ii,23)') 
                self.write_line('fcc(aux,1)=elem_infos(ii,24)')
                
        
                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR')
                self.write_line('*else')    
                self.write_line('*endif')  
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_loc_x_glob_x + '(1) , \',\' , ' + name_loc_x_glob_y + '(1) , \',\' ,' + name_loc_x_glob_z + '(1) , \',\' ,'+ name_loc_y_glob_x + '(1) , \',\' ,' + name_loc_y_glob_y + '(1) , \',\' ,'+ name_loc_y_glob_z + '(1) , \',\' ,' + name_elem_typ + '(1) ')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname_1 + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_ecu + '(1) , \',\' , ' + name_kE + '(1) , \',\' ,' + name_kriss + '(1) , \',\' ,' + name_fcc + '(1) ')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('!')
                self.write_line('!')
                cFile.close()                   


                # Top Layer (nn)
                # --------------------------------------------
                fname = str(step_name) + '_' + 'eps_x_y_xy_top'
                name_ = 'eps_x_y_xy_GP'
                name_elem_nr = 'elem_nr'
                name_eps_x = 'eps_x'
                name_eps_y = 'eps_y'
                name_eps_xy = 'eps_xy'  
                name_bruch = 'eps_bruch'              
                name_coor_intp_layer_x ='coor_intp_layer_x'
                name_coor_intp_layer_y ='coor_intp_layer_y'
                name_coor_intp_layer_z ='coor_intp_layer_z'

        
                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write strains in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,eps_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_xy,ARRAY,4*NrE,1')                
                self.write_line('*DIM,eps_bruch,ARRAY,4*NrE,1') 

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                self.write_line('*DIM,coor_intp_layer_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_z,ARRAY,4*NrE,1')               
  
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')          
                
                # (NrT: Shell ID, NrL=Numbers of Layer) 
                self.write_line('*GET, NrT, ELEM,N_E(ii), attr, secn') # gibt zu einem Element zugehorgie secnum
                self.write_line('*GET, NrL, SHEL, NrT, Prop,NLAY ') # NrT is equal to secnum (not Element number)
                            

                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')                           

                self.write_line('LAYER,NrL')
                self.write_line('elem_nr(aux,1)=N_E(ii)')
                self.write_line('*GET,usedmodel_check,NODE,N_N(kk),SVAR,1 ')
                                              
                self.write_line('*GET,eps_x(aux,1),NODE,N_N(kk),SVAR,69')                
                self.write_line('*GET,eps_y(aux,1),NODE,N_N(kk),SVAR,70')                
                self.write_line('*GET,eps_xy(aux,1),NODE,N_N(kk),SVAR,71')                 
                self.write_line('*GET,eps_bruch(aux,1),NODE,N_N(kk),SVAR,2') 

                self.write_line('*GET,coor_intp_layer_x(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_intp_layer_y(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_intp_layer_z(aux,1),NODE,N_N(kk),SVAR,65')               

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR')
                self.write_line('*else')    
                self.write_line('*endif')  
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')


                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_eps_x + '(1) , \',\' ,' + name_eps_y + '(1) , \',\' ,'+ name_eps_xy + '(1) , \',\' ,' + name_bruch + '(1) , \',\' ,'+ name_coor_intp_layer_x + '(1) , \',\' ,'  + name_coor_intp_layer_y + '(1) , \',\' ,' + name_coor_intp_layer_z + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')
                

                self.write_line('!')
                self.write_line('!')
                cFile.close()       

                
                # Bot Layer (1)
                # --------------------------------------------
                fname = str(step_name) + '_' + 'eps_x_y_xy_bot'
                name_ = 'eps_x_y_xy_GP'
                name_elem_nr = 'elem_nr'
                name_eps_x = 'eps_x'
                name_eps_y = 'eps_y'
                name_eps_xy = 'eps_xy'                
                name_bruch = 'eps_bruch' 
                name_coor_intp_layer_x ='coor_intp_layer_x'
                name_coor_intp_layer_y ='coor_intp_layer_y'
                name_coor_intp_layer_z ='coor_intp_layer_z'

        
                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write strains in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,eps_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_xy,ARRAY,4*NrE,1')  
                self.write_line('*DIM,eps_bruch,ARRAY,4*NrE,1')               

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                self.write_line('*DIM,coor_intp_layer_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_z,ARRAY,4*NrE,1')               
  
                
                # Extract Principal strains from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')                                     

                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')                           

                self.write_line('LAYER,1')
                self.write_line('elem_nr(aux,1)=N_E(ii)')
                self.write_line('*GET,usedmodel_check,NODE,N_N(kk),SVAR,1 ')
                                              
                self.write_line('*GET,eps_x(aux,1),NODE,N_N(kk),SVAR,69')                
                self.write_line('*GET,eps_y(aux,1),NODE,N_N(kk),SVAR,70')                
                self.write_line('*GET,eps_xy(aux,1),NODE,N_N(kk),SVAR,71')    
                self.write_line('*GET,eps_bruch(aux,1),NODE,N_N(kk),SVAR,2')              


                self.write_line('*GET,coor_intp_layer_x(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_intp_layer_y(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_intp_layer_z(aux,1),NODE,N_N(kk),SVAR,65')               

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR')
                self.write_line('*else')    
                self.write_line('*endif')  
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')

                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_eps_x + '(1) , \',\' ,' + name_eps_y + '(1) , \',\' ,'+ name_eps_xy + '(1) , \',\' ,' + name_bruch + '(1) , \',\' ,'+ name_coor_intp_layer_x + '(1) , \',\' ,'  + name_coor_intp_layer_y + '(1) , \',\' ,' + name_coor_intp_layer_z + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('!')
                self.write_line('!')
                cFile.close()                 
          
            else:
                pass 

            # Write principal strains at top and bottom Elementlayer at every GP
            # ------------------------------------------------------------------               
            if 'eps' in fields or 'all' in fields:

                # Top Layer (nn)
                # --------------------------------------------
                fname = str(step_name) + '_' + 'strains_top'
                name_ = 'eps_GP'
                name_elem_nr = 'elem_nr'
                name_eps_1 = 'eps_1'
                name_eps_3 = 'eps_3'
                name_coor_intp_layer_x ='coor_intp_layer_x'
                name_coor_intp_layer_y ='coor_intp_layer_y'
                name_coor_intp_layer_z ='coor_intp_layer_z'

        
                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write strains in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,eps_1,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_3,ARRAY,4*NrE,1')

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                self.write_line('*DIM,coor_intp_layer_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_z,ARRAY,4*NrE,1')               
  
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')          
                
                # (NrT: Shell ID, NrL=Numbers of Layer) 
                self.write_line('*GET, NrT, ELEM,N_E(ii), attr, secn') # gibt zu einem Element zugehorgie secnum
                self.write_line('*GET, NrL, SHEL, NrT, Prop,NLAY ') # NrT is equal to secnum (not Element number)
                            

                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')

                self.write_line('LAYER,NrL')
                self.write_line('elem_nr(aux,1)=N_E(ii)')                                
                self.write_line('*GET,eps_1(aux,1),NODE,N_N(kk),SVAR,5')
                self.write_line('*GET,eps_3(aux,1),NODE,N_N(kk),SVAR,6')
                self.write_line('*GET,coor_intp_layer_x(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_intp_layer_y(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_intp_layer_z(aux,1),NODE,N_N(kk),SVAR,65')               

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR') 
                self.write_line('*else')    
                self.write_line('*endif')                      
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')

                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_eps_1 + '(1) , \',\' ,' + name_eps_3 + '(1) , \',\' ,'+ name_coor_intp_layer_x + '(1) , \',\' ,'  + name_coor_intp_layer_y + '(1) , \',\' ,' + name_coor_intp_layer_z + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('!')
                self.write_line('!')
                cFile.close()      

                # Bot Layer (1)
                # --------------------------------------------
                
                fname = str(step_name) + '_' + 'strains_bot'
                name_ = 'eps_GP'
                name_elem_nr = 'elem_nr'
                name_eps_1 = 'eps_1'
                name_eps_3 = 'eps_3'
                name_coor_intp_layer_x ='coor_intp_layer_x'
                name_coor_intp_layer_y ='coor_intp_layer_y'
                name_coor_intp_layer_z ='coor_intp_layer_z'

        
                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write strains in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,eps_1,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_3,ARRAY,4*NrE,1')

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                self.write_line('*DIM,coor_intp_layer_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_z,ARRAY,4*NrE,1')               
  
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')          
                
                # (NrT: Shell ID, NrL=Numbers of Layer) 
                self.write_line('*GET, NrT, ELEM,N_E(ii), attr, secn') # gibt zu einem Element zugehorgie secnum
                self.write_line('*GET, NrL, SHEL, NrT, Prop,NLAY ') # NrT is equal to secnum (not Element number)
                            

                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')

                self.write_line('LAYER,1')
                self.write_line('elem_nr(aux,1)=N_E(ii)')                                
                self.write_line('*GET,eps_1(aux,1),NODE,N_N(kk),SVAR,5')
                self.write_line('*GET,eps_3(aux,1),NODE,N_N(kk),SVAR,6')
                self.write_line('*GET,coor_intp_layer_x(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_intp_layer_y(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_intp_layer_z(aux,1),NODE,N_N(kk),SVAR,65')               

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR') 
                self.write_line('*else')    
                self.write_line('*endif')                      
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')

                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_eps_1 + '(1) , \',\' ,' + name_eps_3 + '(1) , \',\' ,'+ name_coor_intp_layer_x + '(1) , \',\' ,'  + name_coor_intp_layer_y + '(1) , \',\' ,' + name_coor_intp_layer_z + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('!')
                self.write_line('!')
                cFile.close()    
          
            else:
                pass 


            #  Write Eps 06d top and bot for shear verification
            # ------------------------------------------------------------------
            if 'eps' in fields or 'all' in fields:

                # eps x,y, xy 06 bot 
                # --------------------------------------------
                fname = str(step_name) + '_' + 'eps_x_y_xy_06_bot'
                name_ = 'eps_x_y_xy_06_bot_GP'
                name_elem_nr = 'elem_nr'
                name_eps_x = 'eps_x'
                name_eps_y = 'eps_y'
                name_eps_xy = 'eps_xy'                            
                name_coor_intp_layer_x ='coor_intp_layer_x'
                name_coor_intp_layer_y ='coor_intp_layer_y'
                name_coor_intp_layer_z ='coor_intp_layer_z'
                name_layer_06d_bot ='layer_06d_bot'
        
                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write strains at 06 d from bot in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,eps_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_xy,ARRAY,4*NrE,1')                                

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                self.write_line('*DIM,coor_intp_layer_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_z,ARRAY,4*NrE,1')               
                self.write_line('*DIM,layer_06d_bot,ARRAY,4*NrE,1')       
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')          
                
                # actuel Layer number to exctract eps
                self.write_line('LNr_d06_bot_act=elem_infos(ii,31)')                

                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')                           

                self.write_line('LAYER,LNr_d06_bot_act')
                self.write_line('elem_nr(aux,1)=N_E(ii)')                
                self.write_line('layer_06d_bot(aux,1)=LNr_d06_bot_act')

                self.write_line('*GET,eps_x(aux,1),NODE,N_N(kk),SVAR,69')                
                self.write_line('*GET,eps_y(aux,1),NODE,N_N(kk),SVAR,70')                
                self.write_line('*GET,eps_xy(aux,1),NODE,N_N(kk),SVAR,71')                                                                    

                self.write_line('*GET,coor_intp_layer_x(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_intp_layer_y(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_intp_layer_z(aux,1),NODE,N_N(kk),SVAR,65')               

                self.write_line('*ENDDO')
                
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR')
                self.write_line('*else')    
                self.write_line('*endif')  
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')


                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_eps_x + '(1) , \',\' ,' + name_eps_y + '(1) , \',\' ,'+ name_eps_xy + '(1) , \',\' ,' + name_layer_06d_bot + '(1) , \',\' ,'+ name_coor_intp_layer_x + '(1) , \',\' ,'  + name_coor_intp_layer_y + '(1) , \',\' ,' + name_coor_intp_layer_z + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')
                

                self.write_line('!')
                self.write_line('!')
                cFile.close()       

                # eps x,y, xy 06 top 
                # --------------------------------------------
                fname = str(step_name) + '_' + 'eps_x_y_xy_06_top'
                name_ = 'eps_x_y_xy_06_top_GP'
                name_elem_nr = 'elem_nr'
                name_eps_x = 'eps_x'
                name_eps_y = 'eps_y'
                name_eps_xy = 'eps_xy'                            
                name_coor_intp_layer_x ='coor_intp_layer_x'
                name_coor_intp_layer_y ='coor_intp_layer_y'
                name_coor_intp_layer_z ='coor_intp_layer_z'
                name_layer_06d_top ='layer_06d_top'

        
                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write strains at 06 d from top in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,eps_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,eps_xy,ARRAY,4*NrE,1')                                

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                self.write_line('*DIM,coor_intp_layer_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_intp_layer_z,ARRAY,4*NrE,1')               
                self.write_line('*DIM,layer_06d_top,ARRAY,4*NrE,1')       
  
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')          
                
                # actuel Layer number to exctract eps
                self.write_line('LNr_d06_top_act=elem_infos(ii,32)')                

                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')                           

                self.write_line('LAYER,LNr_d06_top_act')
                self.write_line('elem_nr(aux,1)=N_E(ii)')                
                self.write_line('layer_06d_top(aux,1)=LNr_d06_top_act')  
                                              
                self.write_line('*GET,eps_x(aux,1),NODE,N_N(kk),SVAR,69')                
                self.write_line('*GET,eps_y(aux,1),NODE,N_N(kk),SVAR,70')                
                self.write_line('*GET,eps_xy(aux,1),NODE,N_N(kk),SVAR,71')                                 

                self.write_line('*GET,coor_intp_layer_x(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_intp_layer_y(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_intp_layer_z(aux,1),NODE,N_N(kk),SVAR,65')               

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR')
                self.write_line('*else')    
                self.write_line('*endif')  
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')


                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_eps_x + '(1) , \',\' ,' + name_eps_y + '(1) , \',\' ,'+ name_eps_xy + '(1) , \',\' ,'+ name_layer_06d_top + '(1) , \',\' ,'  + name_coor_intp_layer_x + '(1) , \',\' ,'  + name_coor_intp_layer_y + '(1) , \',\' ,' + name_coor_intp_layer_z + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')
                

                self.write_line('!')
                self.write_line('!')
                cFile.close() 

            # Write steel stresses at reinforcement Layer 1-4 at every GP
            # ------------------------------------------------------------------               
            if 'sig_sr' in fields or 'all' in fields:



                # Benotigte Element infos abspeichern fur steel stresses (MAYBE TO REMOVE)
                # --------------------------------------------
                fname = str(step_name) + '_' + 'sig_sr_elem_infos'
                name_ = 'sig_sr_info_GP'
                name_loc_x_glob_x ='loc_x_glob_x'
                name_loc_x_glob_y ='loc_x_glob_y'
                name_loc_x_glob_z ='loc_x_glob_z'
                name_loc_y_glob_x ='loc_y_glob_x'
                name_loc_y_glob_y ='loc_y_glob_y'
                name_loc_y_glob_z ='loc_y_glob_z'
                name_elem_typ ='elem_typ'

                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write element infos for steel stresses in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,loc_x_glob_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_x_glob_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_x_glob_z,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_y_glob_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_y_glob_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,loc_y_glob_z,ARRAY,4*NrE,1')
                self.write_line('*DIM,elem_typ,ARRAY,4*NrE,1')

                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')                             
                            
                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')


                self.write_line('loc_x_glob_x(aux,1)=elem_infos(ii,3)') 
                self.write_line('loc_x_glob_y(aux,1)=elem_infos(ii,4)') 
                self.write_line('loc_x_glob_z(aux,1)=elem_infos(ii,5)') 
                self.write_line('loc_y_glob_x(aux,1)=elem_infos(ii,6)') 
                self.write_line('loc_y_glob_y(aux,1)=elem_infos(ii,7)') 
                self.write_line('loc_y_glob_z(aux,1)=elem_infos(ii,8)') 
                self.write_line('elem_typ(aux,1)=elem_infos(ii,1)') 
           

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR')
                self.write_line('*else')    
                self.write_line('*endif')  
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')

                self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_loc_x_glob_x + '(1) , \',\' , ' + name_loc_x_glob_y + '(1) , \',\' ,' + name_loc_x_glob_z + '(1) , \',\' ,'+ name_loc_y_glob_x + '(1) , \',\' ,' + name_loc_y_glob_y + '(1) , \',\' ,'+ name_loc_y_glob_z + '(1) , \',\' ,' + name_elem_typ + '(1) ')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('!')
                self.write_line('!')
                cFile.close()                   

                # Reinforcement Layer 1
                # --------------------------------------------
                fname_1L = str(step_name) + '_' + 'sig_sr_1L'
                fname_2L = str(step_name) + '_' + 'sig_sr_2L'
                fname_3L = str(step_name) + '_' + 'sig_sr_3L'
                fname_4L = str(step_name) + '_' + 'sig_sr_4L'
                name_ = 'sig_sr_GP'
                name_elem_nr = 'elem_nr'
                name_sig_sr_1L = 'sig_sr_1L'                 
                name_sig_sr_2L = 'sig_sr_2L'                 
                name_sig_sr_3L = 'sig_sr_3L'                 
                name_sig_sr_4L = 'sig_sr_4L'                 
                
                name_coor_x_sig_sr_1L ='coor_x_sig_sr_1L'
                name_coor_y_sig_sr_1L ='coor_y_sig_sr_1L'
                name_coor_z_sig_sr_1L ='coor_z_sig_sr_1L'
                name_coor_x_sig_sr_2L ='coor_x_sig_sr_2L'
                name_coor_y_sig_sr_2L ='coor_y_sig_sr_2L'
                name_coor_z_sig_sr_2L ='coor_z_sig_sr_2L'
                name_coor_x_sig_sr_3L ='coor_x_sig_sr_3L'
                name_coor_y_sig_sr_3L ='coor_y_sig_sr_3L'
                name_coor_z_sig_sr_3L ='coor_z_sig_sr_3L'    
                name_coor_x_sig_sr_4L ='coor_x_sig_sr_4L'
                name_coor_y_sig_sr_4L ='coor_y_sig_sr_4L'
                name_coor_z_sig_sr_4L ='coor_z_sig_sr_4L'     
                name_psi_1L ='psi_1L'
                name_psi_2L ='psi_2L'
                name_psi_3L ='psi_3L'
                name_psi_4L ='psi_4L'
                name_fsy_1L ='fsy_1L'
                name_fsy_2L ='fsy_2L'
                name_fsy_3L ='fsy_3L'
                name_fsy_4L ='fsy_4L'
                name_fsu_1L ='fsu_1L'
                name_fsu_2L ='fsu_2L'
                name_fsu_3L ='fsu_3L'
                name_fsu_4L ='fsu_4L'

        
                cFile = open(os.path.join(path, filename), 'a')
                self.blank_line()
                self.write_line('! Write steel stresses at crackes in Shell Elements')
                self.blank_line()
                
                # Liste mit allen Elementen aufbauen                 
                self.write_line('allsel')
                self.write_line('nsel,all')
                self.write_line('*get,NrE,elem,0,count') # NrE=Anzahl Elemente
                self.write_line('*dim,N_E,array,NrE,1')
                self.write_line('*vget,N_E,elem,,elist') # N_E=Element liste

                # Aufbau der Arrays
                self.write_line('*DIM,sig_sr_1L,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_2L,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_3L,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_4L,ARRAY,4*NrE,1')
    
                self.write_line('*DIM,elem_nr,ARRAY,4*NrE,1')
                self.write_line('*dim,' + name_ + ', ,4*NrE')
                self.write_line('*DIM,coor_x_sig_sr_1L,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_y_sig_sr_1L,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_z_sig_sr_1L,ARRAY,4*NrE,1')   
                self.write_line('*DIM,coor_x_sig_sr_2L,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_y_sig_sr_2L,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_z_sig_sr_2L,ARRAY,4*NrE,1')    
                self.write_line('*DIM,coor_x_sig_sr_3L,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_y_sig_sr_3L,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_z_sig_sr_3L,ARRAY,4*NrE,1') 
                self.write_line('*DIM,coor_x_sig_sr_4L,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_y_sig_sr_4L,ARRAY,4*NrE,1')
                self.write_line('*DIM,coor_z_sig_sr_4L,ARRAY,4*NrE,1')     
                self.write_line('*DIM,psi_1L,ARRAY,4*NrE,1')
                self.write_line('*DIM,psi_2L,ARRAY,4*NrE,1')
                self.write_line('*DIM,psi_3L,ARRAY,4*NrE,1')
                self.write_line('*DIM,psi_4L,ARRAY,4*NrE,1')                                                 
                self.write_line('*DIM,fsy_1L,ARRAY,4*NrE,1')
                self.write_line('*DIM,fsy_2L,ARRAY,4*NrE,1')
                self.write_line('*DIM,fsy_3L,ARRAY,4*NrE,1')
                self.write_line('*DIM,fsy_4L,ARRAY,4*NrE,1')  
                self.write_line('*DIM,fsu_1L,ARRAY,4*NrE,1')
                self.write_line('*DIM,fsu_2L,ARRAY,4*NrE,1')
                self.write_line('*DIM,fsu_3L,ARRAY,4*NrE,1')
                self.write_line('*DIM,fsu_4L,ARRAY,4*NrE,1')                  
                
                # Extract Principal stresses from usermat
                self.write_line('aux=0')

                self.write_line('*DO,ii,1,NrE') # Loop uber alle Elemente                
                self.write_line('ESEL,S,ELEM, ,N_E(ii)')
                
                self.write_line('*if,elem_infos(ii,1),EQ,1,THEN  ')          
                                           

                self.write_line('NSLE,ALL')
                self.write_line('*GET,NrN,NODE,0,COUNT ')
                self.write_line('*DIM,N_N,ARRAY,NrN,1')
                self.write_line('*VGET,N_N,NODE, ,NLIST')
                self.write_line('*DO,kk,1,NrN')
                self.write_line('aux = aux+1')

                self.write_line('*GET,layer_nr_1L,NODE,N_N(kk),SVAR,58')
                self.write_line('*GET,layer_nr_2L,NODE,N_N(kk),SVAR,59')
                self.write_line('*GET,layer_nr_3L,NODE,N_N(kk),SVAR,60')
                self.write_line('*GET,layer_nr_4L,NODE,N_N(kk),SVAR,61')
                self.write_line('elem_nr(aux,1)=N_E(ii)') 

                # Layer 1 for x- and y direction
                self.write_line('LAYER,layer_nr_1L')
                self.write_line('*GET,sig_sr_1L_Druck_x,NODE,N_N(kk),SVAR,13')
                self.write_line('*GET,sig_sr_1L_Druckfeld_x,NODE,N_N(kk),SVAR,20')
                self.write_line('*GET,sig_sr_1L_DruckZug_x,NODE,N_N(kk),SVAR,34')
                self.write_line('*GET,sig_sr_1L_Zug_x,NODE,N_N(kk),SVAR,40')
                self.write_line('sig_sr_1L_x=sig_sr_1L_Druck_x+sig_sr_1L_Druckfeld_x+sig_sr_1L_DruckZug_x+sig_sr_1L_Zug_x')
                self.write_line('*GET,sig_sr_1L_Druck_y,NODE,N_N(kk),SVAR,14')
                self.write_line('*GET,sig_sr_1L_Druckfeld_y,NODE,N_N(kk),SVAR,21')
                self.write_line('*GET,sig_sr_1L_DruckZug_y,NODE,N_N(kk),SVAR,35')
                self.write_line('*GET,sig_sr_1L_Zug_y,NODE,N_N(kk),SVAR,41')
                self.write_line('sig_sr_1L_y=sig_sr_1L_Druck_y+sig_sr_1L_Druckfeld_y+sig_sr_1L_DruckZug_y+sig_sr_1L_Zug_y')
                self.write_line('sig_sr_1L(aux,1)=sig_sr_1L_x+sig_sr_1L_y') # fiktive x und y Richtung zusammenzahlen
                self.write_line('*GET,coor_x_sig_sr_1L(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_y_sig_sr_1L(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_z_sig_sr_1L(aux,1),NODE,N_N(kk),SVAR,65')                  
                self.write_line('psi_1L(aux,1)=elem_infos(ii,9)')           
                self.write_line('fsy_1L(aux,1)=elem_infos(ii,13)') 
                self.write_line('fsu_1L(aux,1)=elem_infos(ii,17)') 

                # Layer 2 for x- and y direction
                self.write_line('LAYER,layer_nr_2L')
                self.write_line('*GET,sig_sr_2L_Druck_x,NODE,N_N(kk),SVAR,13')
                self.write_line('*GET,sig_sr_2L_Druckfeld_x,NODE,N_N(kk),SVAR,20')
                self.write_line('*GET,sig_sr_2L_DruckZug_x,NODE,N_N(kk),SVAR,34')
                self.write_line('*GET,sig_sr_2L_Zug_x,NODE,N_N(kk),SVAR,40')
                self.write_line('sig_sr_2L_x=sig_sr_2L_Druck_x+sig_sr_2L_Druckfeld_x+sig_sr_2L_DruckZug_x+sig_sr_2L_Zug_x')
                self.write_line('*GET,sig_sr_2L_Druck_y,NODE,N_N(kk),SVAR,14')
                self.write_line('*GET,sig_sr_2L_Druckfeld_y,NODE,N_N(kk),SVAR,21')
                self.write_line('*GET,sig_sr_2L_DruckZug_y,NODE,N_N(kk),SVAR,35')
                self.write_line('*GET,sig_sr_2L_Zug_y,NODE,N_N(kk),SVAR,41')
                self.write_line('sig_sr_2L_y=sig_sr_2L_Druck_y+sig_sr_2L_Druckfeld_y+sig_sr_2L_DruckZug_y+sig_sr_2L_Zug_y')
                self.write_line('sig_sr_2L(aux,1)=sig_sr_2L_x+sig_sr_2L_y') # fiktive x und y Richtung zusammenzahlen
                self.write_line('*GET,coor_x_sig_sr_2L(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_y_sig_sr_2L(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_z_sig_sr_2L(aux,1),NODE,N_N(kk),SVAR,65')   
                self.write_line('psi_2L(aux,1)=elem_infos(ii,10)')  
                self.write_line('fsy_2L(aux,1)=elem_infos(ii,14)')  
                self.write_line('fsu_2L(aux,1)=elem_infos(ii,18)') 

                # Layer 3 for x- and y direction
                self.write_line('LAYER,layer_nr_3L')
                self.write_line('*GET,sig_sr_3L_Druck_x,NODE,N_N(kk),SVAR,13')
                self.write_line('*GET,sig_sr_3L_Druckfeld_x,NODE,N_N(kk),SVAR,20')
                self.write_line('*GET,sig_sr_3L_DruckZug_x,NODE,N_N(kk),SVAR,34')
                self.write_line('*GET,sig_sr_3L_Zug_x,NODE,N_N(kk),SVAR,40')
                self.write_line('sig_sr_3L_x=sig_sr_3L_Druck_x+sig_sr_3L_Druckfeld_x+sig_sr_3L_DruckZug_x+sig_sr_3L_Zug_x')
                self.write_line('*GET,sig_sr_3L_Druck_y,NODE,N_N(kk),SVAR,14')
                self.write_line('*GET,sig_sr_3L_Druckfeld_y,NODE,N_N(kk),SVAR,21')
                self.write_line('*GET,sig_sr_3L_DruckZug_y,NODE,N_N(kk),SVAR,35')
                self.write_line('*GET,sig_sr_3L_Zug_y,NODE,N_N(kk),SVAR,41')
                self.write_line('sig_sr_3L_y=sig_sr_3L_Druck_y+sig_sr_3L_Druckfeld_y+sig_sr_3L_DruckZug_y+sig_sr_3L_Zug_y')
                self.write_line('sig_sr_3L(aux,1)=sig_sr_3L_x+sig_sr_3L_y') # fiktive x und y Richtung zusammenzahlen
                self.write_line('*GET,coor_x_sig_sr_3L(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_y_sig_sr_3L(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_z_sig_sr_3L(aux,1),NODE,N_N(kk),SVAR,65')     
                self.write_line('psi_3L(aux,1)=elem_infos(ii,11)')   
                self.write_line('fsy_3L(aux,1)=elem_infos(ii,15)') 
                self.write_line('fsu_3L(aux,1)=elem_infos(ii,19)')  

                # Layer 4 for x- and y direction
                self.write_line('LAYER,layer_nr_4L')
                self.write_line('*GET,sig_sr_4L_Druck_x,NODE,N_N(kk),SVAR,13')
                self.write_line('*GET,sig_sr_4L_Druckfeld_x,NODE,N_N(kk),SVAR,20')
                self.write_line('*GET,sig_sr_4L_DruckZug_x,NODE,N_N(kk),SVAR,34')
                self.write_line('*GET,sig_sr_4L_Zug_x,NODE,N_N(kk),SVAR,40')
                self.write_line('sig_sr_4L_x=sig_sr_4L_Druck_x+sig_sr_4L_Druckfeld_x+sig_sr_4L_DruckZug_x+sig_sr_4L_Zug_x')
                self.write_line('*GET,sig_sr_4L_Druck_y,NODE,N_N(kk),SVAR,14')
                self.write_line('*GET,sig_sr_4L_Druckfeld_y,NODE,N_N(kk),SVAR,21')
                self.write_line('*GET,sig_sr_4L_DruckZug_y,NODE,N_N(kk),SVAR,35')
                self.write_line('*GET,sig_sr_4L_Zug_y,NODE,N_N(kk),SVAR,41')
                self.write_line('sig_sr_4L_y=sig_sr_4L_Druck_y+sig_sr_4L_Druckfeld_y+sig_sr_4L_DruckZug_y+sig_sr_4L_Zug_y')
                self.write_line('sig_sr_4L(aux,1)=sig_sr_4L_x+sig_sr_4L_y') # fiktive x und y Richtung zusammenzahlen
                self.write_line('*GET,coor_x_sig_sr_4L(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_y_sig_sr_4L(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_z_sig_sr_4L(aux,1),NODE,N_N(kk),SVAR,65') 
                self.write_line('psi_4L(aux,1)=elem_infos(ii,12)')    
                self.write_line('fsy_4L(aux,1)=elem_infos(ii,16)')  
                self.write_line('fsu_4L(aux,1)=elem_infos(ii,20)')         

                self.write_line('*ENDDO')
                self.write_line('*DEL,N_N,,NOPR')
                self.write_line('*DEL,NrN,,NOPR') 
                self.write_line('*else')    
                self.write_line('*endif')                      
                self.write_line('*ENDDO')
                self.write_line('*DEL,NrE,,NOPR')
                self.write_line('*DEL,N_E,,NOPR')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname_1L + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_sr_1L + '(1) , \',\' , ' + name_coor_x_sig_sr_1L + '(1) , \',\' ,'+ name_coor_y_sig_sr_1L + '(1) , \',\' ,'  + name_coor_z_sig_sr_1L + '(1), \',\' ,'  + name_psi_1L + '(1), \',\' ,'  + name_fsy_1L + '(1), \',\' ,'  + name_fsu_1L + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname_2L + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_sr_2L + '(1) , \',\' ,' + name_coor_x_sig_sr_2L + '(1) , \',\' ,'+ name_coor_y_sig_sr_2L + '(1) , \',\' ,'  + name_coor_z_sig_sr_2L + '(1), \',\' ,'  + name_psi_2L + '(1), \',\' ,'  + name_fsy_2L + '(1), \',\' ,'  + name_fsu_2L + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')                

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname_3L + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_sr_3L + '(1) , \',\' , ' + name_coor_x_sig_sr_3L + '(1) , \',\' ,'+ name_coor_y_sig_sr_3L + '(1) , \',\' ,'  + name_coor_z_sig_sr_3L + '(1), \',\' ,'  + name_psi_3L + '(1), \',\' ,'  + name_fsy_3L + '(1), \',\' ,'  + name_fsu_3L + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')   

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname_4L + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_sr_4L + '(1) , \',\' , '  + name_coor_x_sig_sr_4L + '(1) , \',\' ,'+ name_coor_y_sig_sr_4L + '(1) , \',\' ,'  + name_coor_z_sig_sr_4L + '(1), \',\' ,'  + name_psi_4L + '(1), \',\' ,'  + name_fsy_4L + '(1), \',\' ,'  + name_fsu_4L + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')   

                self.write_line('!')
                self.write_line('!')
                cFile.close()      

            else:
                pass             

        self.blank_line()
        self.blank_line()
        self.write_line('*cfopen,run_ansys_check,txt ')
        self.blank_line()
        self.blank_line()
    