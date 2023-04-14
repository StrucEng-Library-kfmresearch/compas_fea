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
            
            # Write Displacement at nodes
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

            # Write Shell forces and moments (only for Shell elements)
            if 'sf' in fields or 'all' in fields:

                fname = str(step_name) + '_' + 'shell_forces_moments'
                name_ = 'ele_sf'


                elements = self.structure.elements
                num_elements = len(elements)
                properties = self.structure.element_properties
                sections = self.structure.sections
                sets = self.structure.sets

                self.write_line('*DIM,elem_infos,ARRAY,{0},2'.format(num_elements))   #set kp origin 

                # Bestimmung des zur einem Element (number) gehoerende Elementtyp (shelle or others)
                for key in sorted(properties):
                    property = properties[key]
                    section = sections[property.section]
                    stype = section.__name__
                    elset = property.elset
                    selection = property.elements if property.elements else sets[elset].selection # Zum elemeset zugehorige Elementnummern
    
                    if stype == 'ShellSection':
                        stype_number=1
                    else:
                        stype_number=0
                    
                    for ele_num in sorted(selection):                    
                        ele_num_adj=ele_num+1
                        self.write_line('elem_infos({0},1)={1} ! Elementyp e.g. shell=0, mpc and other =1'.format(ele_num_adj,stype_number))   # Elementype     
                        self.write_line('elem_infos({0},2)={1} ! Elementnumber'.format(ele_num_adj,ele_num_adj))   # Elementnumber 

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
                self.write_line('*if,elem_infos(i,1),EQ,1,THEN  ')             
                self.write_line('*get, eforces(i, 1), ETAB, 1, ELEM, enum(i) ') # N11 In-plane forces (per unit length)
                self.write_line('*get, eforces(i, 2), ETAB, 2, ELEM, enum(i) ') # N22 In-plane forces (per unit length)
                self.write_line('*get, eforces(i, 3), ETAB, 3, ELEM, enum(i) ') # N12 In-plane forces (per unit length)
                self.write_line('*get, eforces(i, 4), ETAB, 7, ELEM, enum(i) ') # Q13 Out-of-plane moments (per unit length)
                self.write_line('*get, eforces(i, 5), ETAB, 8, ELEM, enum(i) ') # Q23 Out-of-plane moments (per unit length)
                self.write_line('*get, eforces(i, 6), ETAB, 4, ELEM, enum(i) ') # M11 Out-of-plane moments (per unit length)
                self.write_line('*get, eforces(i, 7), ETAB, 5, ELEM, enum(i) ') # M22 Out-of-plane moments (per unit length)
                self.write_line('*get, eforces(i, 8), ETAB, 6, ELEM, enum(i) ') # M12 Out-of-plane moments (per unit length)       
                self.write_line('eforces(i, 9)=elem_infos(i,1)') # M12 Out-of-plane moments (per unit length)    
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
        
                #self.write_line('*cfopen,' + out_path + '/' + fname + ',txt \n')
                #print('slayer')
                #self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                #self.write_line('*CFWRITE, stresses, Enumber, sf1, sf2, sf3 \n')
                #self.write_line('*do,i,1,nelem \n')
                #self.write_line('*CFWRITE, stresses, enum(i,1), eforces(i,1), eforces(i,2), eforces(i,3) \n')
                #self.write_line('*vwrite, ' + name_ + '(1) , \',\'  , ' + name_x + '(1) , \',\' , ' + name_y + '(1) , \',\' ,' + name_z + '(1)')
                #self.write_line('*CFWRITE, ' +name_ + '(1) , ' + 'eforces(i,1), eforces(i,2), eforces(i,3) \n')
                #self.write_line('*Enddo \n')
                #print('SCHLACH')
                #self.write_line('ESEL, ALL \n')
                #self.write_line('ETABLE, ERAS \n')
                #self.write_line('! \n')
                #cFile.close()
                #self.write_line('*vfill,' + name_ + '(1),ramp,1,1')
                #self.write_line('*cfopen,' + out_path + '/' + fname + ',txt')
                #self.write_line('*do,i,1,nelem \n')            
                #self.write_line('*vwrite, ' + 'eforces' + '(1,i) , \',\' , ' + 'eforces' + '(1,i) , \',\' ,' + 'eforces' + '(1,i)')
                #self.write_line('(          F4.0,       A,       ES,           A,          ES,          A,      ES)')
                #self.write_line('*Enddo \n')
                #self.write_line('!')
                #self.write_line('!')
                #cFile.close()


                #self.blank_line()
                #self.write_line('! Write Shelll Forces')
                #self.blank_line()
                #self.write_line('allsel')
                #self.write_line('nsel,all')
                #self.write_line('*GET,NrE,ELEM,0,COUNT ')
                #self.write_line('*DIM,N_E,ARRAY,NrE,1')
                #self.write_line('*VGET,N_E,ELEM, ,ELIST')
                #self.write_line('ETAB,ERASE')
                #self.blank_line()            
                #self.write_line('*DIM,elem_nr,ARRAY,NrE,1')
                #self.write_line('*DIM,mvn,ARRAY,NrE,8')
                #self.blank_line()
                #self.write_line('*DO,ii,1,NrE')
                #self.write_line('elem_nr(ii,1) = N_E(ii)')
                #self.write_line('ETABLE,NX,SMISC,1')                   
                #self.write_line('*GET,mvn(ii,1),ETAB,1,ELEM,N_E(ii)')   
                #self.write_line('ETABLE,NY,SMISC,2')                         
                #self.write_line('*GET,mvn(ii,2),ETAB,2,ELEM,N_E(ii)')   
                #self.write_line('ETABLE,NXY,SMISC,3')                         
                #self.write_line('*GET,mvn(ii,3),ETAB,3,ELEM,N_E(ii)')   
                #self.write_line('ETABLE,MX,SMISC,4')                  
                #self.write_line('*GET,mvn(ii,4),ETAB,4,ELEM,N_E(ii)')   
                #self.write_line('ETABLE,MY,SMISC,5')  
                #self.write_line('*GET,mvn(ii,5),ETAB,5,ELEM,N_E(ii)')   
                #self.write_line('ETABLE,MXY,SMISC,6')                   
                #self.write_line('*GET,mvn(ii,6),ETAB,6,ELEM,N_E(ii)')   
                #self.write_line('ETABLE,VX,SMISC,7')  
                #self.write_line('*GET,mvn(ii,7),ETAB,7,ELEM,N_E(ii)')   
                #self.write_line('ETABLE,VY,SMISC,8')  
                #self.write_line('*GET,mvn(ii,8),ETAB,8,ELEM,N_E(ii)')  
                #self.write_line('*ENDDO')
                #self.blank_line()
                
    

            
                #self.blank_line()
                #self.write_line('*DEL,NrN,,NOPR')
                #self.write_line('*DEL,N_N,,NOPR')
                #self.write_line('*DEL,NrE,,NOPR')
                #self.write_line('*DEL,N_E,,NOPR')
                #self.blank_line()


                #self.write_line('Arg1' + '=' + "'"+ 'mvn' + "'")
                #self.write_line('Arg2' + '=' + "'"+ 'f20.5' + "'")
                #self.write_line('Arg3=13')
                #self.write_line('Arg4' + '=' + "'"+ 'eforces.txt' + "'")

                #self.write_line('~eui,' + "'"+  'set arrname [ans_getvalue PARM,Arg1,value]'+ "'")
                #self.write_line('~eui,' + "'"+  'set formatdescriptor [ans_getvalue PARM,Arg2,value]'+ "'")
                #self.write_line('~eui,' + "'"+  'set numcol [ans_getvalue PARM,Arg3,value]'+ "'")
                #self.write_line('~eui,' + "'"+  'set filename [ans_getvalue PARM,Arg4,value]'+ "'")            
                #self.write_line('~eui,' + "'"+  'set arrname [string trim $arrname]'+ "'")
                #self.write_line('~eui,' + "'"+  'set formatdescriptor [string trim $formatdescriptor]'+ "'")
                #self.write_line('~eui,' + "'"+  'set filename [string trim $filename]'+ "'") 
                #self.write_line('~eui,' + "'"+  'set cmd "*mwrite,$arrname"'+ "'")
                #self.write_line('~eui,' + "'"+  'append cmd "(1,1),"'+ "'")
                #self.write_line('~eui,' + "'"+  'append cmd [file rootname $filename]'+ "'")            
                #self.write_line('~eui,' + "'"+  'append cmd ","'+ "'")
                #self.write_line('~eui,' + "'"+  'append cmd [string trimleft [file extension $filename] "."]'+ "'")
                #self.write_line('~eui,' + "'"+  'append cmd' + ' "' + '\\n' + '"' + "'")
                #self.write_line('~eui,' + "'"+  'append cmd "($numcol$formatdescriptor)"' + "'")
                #self.write_line('~eui,' + "'"+  'set fid [open mk_mwritehelper.mac w]'+ "'")
                #self.write_line('~eui,' + "'"+  'puts $fid $cmd'+ "'")
                #self.write_line('~eui,' + "'"+  'close $fid'+ "'")
                #self.write_line('~eui,' + "'"+  'ans_sendcommand mk_mwritehelper'+ "'")
                #self.write_line('~eui,' + "'"+  'file delete mk_mwritehelper.mac'+ "'")

                #self.write_line('*DEL,mvn,,NOPR')

            else:
                pass 

        self.blank_line()
        self.blank_line()
        self.write_line('*cfopen,run_ansys_check,txt ')
        self.blank_line()
        self.blank_line()
            

