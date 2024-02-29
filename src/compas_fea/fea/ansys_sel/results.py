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

        # Write intro part
        self.write_section('Results')
        self.blank_line()
        
        self.write_line('/post1')


        # Schleife uber alle vorhandenen loadsteps

        # Bestimmung der lstep und sbstep welche als txt gespeichert werden sollen
        sbstep=sbstep # Substep
        step_list=structure.steps_order
        
        # Schleife uber alle vorhandenen loadsteps
        for single_lstep in lstep:
            
            lstep_index=step_list.index(single_lstep)         
            
            
            self.write_line('set,{0},{1}'.format(lstep_index,sbstep))
            print('---------------------')
            
            # Define steps for output
            name = structure.name
            path = structure.path
            step_name = structure.steps_order[lstep_index]
            

            # Loschen der Dateien in Output Pfad
            out_path = os.path.join(path, name + '_output')
            dir = out_path
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            else:
                for f in os.listdir(dir):
                    os.remove(os.path.join(dir, f))
    

            filename = name + '_extract.txt'
            
            # General element infos (Werden dann unten im entsprechenden benotigen speicherung der variabeln wiederverwendet)
            # ------------------------------------------------------------------
                      
            elements = self.structure.elements
            num_elements = len(elements)
            properties = self.structure.element_properties
            sections = self.structure.sections
            sets = self.structure.sets
            

            self.write_line('*DIM,elem_infos,ARRAY,{0},8'.format(num_elements))   
            fname = str(step_name) + '_' + 'elem_infos'
            name_ = 'elem_infos'
            
            # Bestimmung der elementspezifischen Werte
            for key in sorted(properties):
                property = properties[key]
                section = sections[property.section]
                stype = section.__name__
                elset = property.elset
                element_set = sets[property.elset]
                selection = property.elements if property.elements else sets[elset].selection # Zum elemeset zugehorige Elementnummern
                
                if stype == 'ShellSection':
                    stype_number=1 # Schalenelement

                else:
                    stype_number=0 # MPC or Other
                    
                for ele_num in sorted(selection):                    
                    ele_num_adj=ele_num+1

                    if stype_number == 1: # Schalenelement
                        loc_coor_EV_XA= section.loc_coords_EV_XA
                        loc_coor_EV_YA= section.loc_coords_EV_YA
                        loc_coor_EV_ZA= section.loc_coords_EV_ZA
                        
                        
                        if element_set.type != 'node':   
                            e_x=loc_coor_EV_XA.get('EV_XA',None)
                            e_y=loc_coor_EV_YA.get('EV_YA',None)
                            e_z=loc_coor_EV_ZA.get('EV_ZA',None)
                        
                        self.write_line('elem_infos({0},1)={1} ! Elementyp e.g. shell=1, mpc and other =0'.format(ele_num_adj,stype_number))   # Elementype     
                        self.write_line('elem_infos({0},2)={1} ! Elementnumber'.format(ele_num_adj,ele_num_adj))   # Elementnumber 
                        self.write_line('elem_infos({0},3)={1}! loccal x-axis direction global x'.format(ele_num_adj,e_x[0]))     
                        self.write_line('elem_infos({0},4)={1}! loccal x-axis direction global y'.format(ele_num_adj,e_x[1]))      
                        self.write_line('elem_infos({0},5)={1}! loccal x-axis direction global z'.format(ele_num_adj,e_x[2]))    
                        self.write_line('elem_infos({0},6)={1}! loccal y-axis direction global x'.format(ele_num_adj,e_y[0]))      
                        self.write_line('elem_infos({0},7)={1}! loccal y-axis direction global y'.format(ele_num_adj,e_y[1])) 
                        self.write_line('elem_infos({0},8)={1}! loccal y-axis direction global z'.format(ele_num_adj,e_y[2]))         
                
                    else:
                        self.write_line('elem_infos({0},1)={1} ! Elementyp e.g. shell=0, mpc and other =1'.format(ele_num_adj,stype_number))   # Elementype     
                        self.write_line('elem_infos({0},2)={1} ! Elementnumber'.format(ele_num_adj,ele_num_adj))   # Elementnumber 
                        self.write_line('elem_infos({0},3)=0'.format(ele_num_adj))
                        self.write_line('elem_infos({0},4)=0'.format(ele_num_adj))
                        self.write_line('elem_infos({0},5)=0'.format(ele_num_adj)) 
                        self.write_line('elem_infos({0},6)=0'.format(ele_num_adj)) 
                        self.write_line('elem_infos({0},7)=0'.format(ele_num_adj)) 
                        self.write_line('elem_infos({0},8)=0'.format(ele_num_adj))
                                       

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

            # Write DISPLACEMENTS at nodes
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

            # Write SHELL FORCES AND MOMENTS (only for Shell elements)
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
        
            # Write STRESSES at top and bottom Elementlayer from usermat results
            # ------------------------------------------------------------------               
            if 's' in fields or 'all' in fields:

                # Benotigte Element infos abspeichern fur stresses
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
                self.write_line('*DIM,testing,ARRAY,4*NrE,1')
  
                
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
                self.write_line('fcc_eff(aux,1)=0')

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
                self.write_line('fcc_eff(aux,1)=0')
                
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
                self.write_line('fcc_eff(aux,1)=0')

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
                self.write_line('fcc_eff(aux,1)=0')
                
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

            # Write STRAINS at top and bottom Elementlayer from usermat results
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

            # Write STEEL STRESSES at reinforcement Layer 1-4 from usermat results
            # ------------------------------------------------------------------               
            if 'sig_sr' in fields or 'all' in fields:

                # Reinforcement Layer 1
                # --------------------------------------------
                fname_1L = str(step_name) + '_' + 'sig_sr_1L'
                fname_2L = str(step_name) + '_' + 'sig_sr_2L'
                fname_3L = str(step_name) + '_' + 'sig_sr_3L'
                fname_4L = str(step_name) + '_' + 'sig_sr_4L'
                name_ = 'sig_sr_GP'
                name_elem_nr = 'elem_nr'
                name_sig_sr_1L_x = 'sig_sr_1L_x'
                name_sig_sr_2L_x = 'sig_sr_2L_x'
                name_sig_sr_3L_x = 'sig_sr_3L_x'
                name_sig_sr_4L_x = 'sig_sr_4L_x'
                name_sig_sr_1L_y = 'sig_sr_1L_y'
                name_sig_sr_2L_y = 'sig_sr_2L_y'
                name_sig_sr_3L_y = 'sig_sr_3L_y'
                name_sig_sr_4L_y = 'sig_sr_4L_y'
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
                self.write_line('*DIM,sig_sr_1L_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_2L_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_3L_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_4L_x,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_1L_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_2L_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_3L_y,ARRAY,4*NrE,1')
                self.write_line('*DIM,sig_sr_4L_y,ARRAY,4*NrE,1')
    
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
                self.write_line('sig_sr_1L_x(aux,1)=sig_sr_1L_Druck_x+sig_sr_1L_Druckfeld_x+sig_sr_1L_DruckZug_x+sig_sr_1L_Zug_x')
                self.write_line('*GET,sig_sr_1L_Druck_y,NODE,N_N(kk),SVAR,14')
                self.write_line('*GET,sig_sr_1L_Druckfeld_y,NODE,N_N(kk),SVAR,21')
                self.write_line('*GET,sig_sr_1L_DruckZug_y,NODE,N_N(kk),SVAR,35')
                self.write_line('*GET,sig_sr_1L_Zug_y,NODE,N_N(kk),SVAR,41')
                self.write_line('sig_sr_1L_y(aux,1)=sig_sr_1L_Druck_y+sig_sr_1L_Druckfeld_y+sig_sr_1L_DruckZug_y+sig_sr_1L_Zug_y')
                self.write_line('*GET,coor_x_sig_sr_1L(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_y_sig_sr_1L(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_z_sig_sr_1L(aux,1),NODE,N_N(kk),SVAR,65')                  
                
                # Layer 2 for x- and y direction
                self.write_line('LAYER,layer_nr_2L')
                self.write_line('*GET,sig_sr_2L_Druck_x,NODE,N_N(kk),SVAR,13')
                self.write_line('*GET,sig_sr_2L_Druckfeld_x,NODE,N_N(kk),SVAR,20')
                self.write_line('*GET,sig_sr_2L_DruckZug_x,NODE,N_N(kk),SVAR,34')
                self.write_line('*GET,sig_sr_2L_Zug_x,NODE,N_N(kk),SVAR,40')
                self.write_line('sig_sr_2L_x(aux,1)=sig_sr_2L_Druck_x+sig_sr_2L_Druckfeld_x+sig_sr_2L_DruckZug_x+sig_sr_2L_Zug_x')
                self.write_line('*GET,sig_sr_2L_Druck_y,NODE,N_N(kk),SVAR,14')
                self.write_line('*GET,sig_sr_2L_Druckfeld_y,NODE,N_N(kk),SVAR,21')
                self.write_line('*GET,sig_sr_2L_DruckZug_y,NODE,N_N(kk),SVAR,35')
                self.write_line('*GET,sig_sr_2L_Zug_y,NODE,N_N(kk),SVAR,41')
                self.write_line('sig_sr_2L_y(aux,1)=sig_sr_2L_Druck_y+sig_sr_2L_Druckfeld_y+sig_sr_2L_DruckZug_y+sig_sr_2L_Zug_y')
                self.write_line('*GET,coor_x_sig_sr_2L(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_y_sig_sr_2L(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_z_sig_sr_2L(aux,1),NODE,N_N(kk),SVAR,65')     

                # Layer 3 for x- and y direction
                self.write_line('LAYER,layer_nr_3L')
                self.write_line('*GET,sig_sr_3L_Druck_x,NODE,N_N(kk),SVAR,13')
                self.write_line('*GET,sig_sr_3L_Druckfeld_x,NODE,N_N(kk),SVAR,20')
                self.write_line('*GET,sig_sr_3L_DruckZug_x,NODE,N_N(kk),SVAR,34')
                self.write_line('*GET,sig_sr_3L_Zug_x,NODE,N_N(kk),SVAR,40')
                self.write_line('sig_sr_3L_x(aux,1)=sig_sr_3L_Druck_x+sig_sr_3L_Druckfeld_x+sig_sr_3L_DruckZug_x+sig_sr_3L_Zug_x')
                self.write_line('*GET,sig_sr_3L_Druck_y,NODE,N_N(kk),SVAR,14')
                self.write_line('*GET,sig_sr_3L_Druckfeld_y,NODE,N_N(kk),SVAR,21')
                self.write_line('*GET,sig_sr_3L_DruckZug_y,NODE,N_N(kk),SVAR,35')
                self.write_line('*GET,sig_sr_3L_Zug_y,NODE,N_N(kk),SVAR,41')
                self.write_line('sig_sr_3L_y(aux,1)=sig_sr_3L_Druck_y+sig_sr_3L_Druckfeld_y+sig_sr_3L_DruckZug_y+sig_sr_3L_Zug_y')
                self.write_line('*GET,coor_x_sig_sr_3L(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_y_sig_sr_3L(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_z_sig_sr_3L(aux,1),NODE,N_N(kk),SVAR,65')     

                # Layer 4 for x- and y direction
                self.write_line('LAYER,layer_nr_4L')
                self.write_line('*GET,sig_sr_4L_Druck_x,NODE,N_N(kk),SVAR,13')
                self.write_line('*GET,sig_sr_4L_Druckfeld_x,NODE,N_N(kk),SVAR,20')
                self.write_line('*GET,sig_sr_4L_DruckZug_x,NODE,N_N(kk),SVAR,34')
                self.write_line('*GET,sig_sr_4L_Zug_x,NODE,N_N(kk),SVAR,40')
                self.write_line('sig_sr_4L_x(aux,1)=sig_sr_4L_Druck_x+sig_sr_4L_Druckfeld_x+sig_sr_4L_DruckZug_x+sig_sr_4L_Zug_x')
                self.write_line('*GET,sig_sr_4L_Druck_y,NODE,N_N(kk),SVAR,14')
                self.write_line('*GET,sig_sr_4L_Druckfeld_y,NODE,N_N(kk),SVAR,21')
                self.write_line('*GET,sig_sr_4L_DruckZug_y,NODE,N_N(kk),SVAR,35')
                self.write_line('*GET,sig_sr_4L_Zug_y,NODE,N_N(kk),SVAR,41')
                self.write_line('sig_sr_4L_y(aux,1)=sig_sr_4L_Druck_y+sig_sr_4L_Druckfeld_y+sig_sr_4L_DruckZug_y+sig_sr_4L_Zug_y')
                self.write_line('*GET,coor_x_sig_sr_4L(aux,1),NODE,N_N(kk),SVAR,63')
                self.write_line('*GET,coor_y_sig_sr_4L(aux,1),NODE,N_N(kk),SVAR,64')
                self.write_line('*GET,coor_z_sig_sr_4L(aux,1),NODE,N_N(kk),SVAR,65')          

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
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_sr_1L_x + '(1) , \',\' , ' + name_sig_sr_1L_y + '(1) , \',\' ,' + name_coor_x_sig_sr_1L + '(1) , \',\' ,'+ name_coor_y_sig_sr_1L + '(1) , \',\' ,'  + name_coor_z_sig_sr_1L + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname_2L + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_sr_2L_x + '(1) , \',\' , ' + name_sig_sr_2L_y + '(1) , \',\' ,' + name_coor_x_sig_sr_2L + '(1) , \',\' ,'+ name_coor_y_sig_sr_2L + '(1) , \',\' ,'  + name_coor_z_sig_sr_2L + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')                

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname_3L + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_sr_3L_x + '(1) , \',\' , ' + name_sig_sr_3L_y + '(1) , \',\' ,' + name_coor_x_sig_sr_3L + '(1) , \',\' ,'+ name_coor_y_sig_sr_3L + '(1) , \',\' ,'  + name_coor_z_sig_sr_3L + '(1)')
                self.write_line('(F100000.0,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES,A,ES)')
                self.write_line('*cfclose \n')   

                self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                self.write_line('*cfopen,' + out_path + '/' + fname_4L + ',txt')
                self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_elem_nr + '(1) , \',\' , ' + name_sig_sr_4L_x + '(1) , \',\' , ' + name_sig_sr_4L_y + '(1) , \',\' ,' + name_coor_x_sig_sr_4L + '(1) , \',\' ,'+ name_coor_y_sig_sr_4L + '(1) , \',\' ,'  + name_coor_z_sig_sr_4L + '(1)')
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
            

