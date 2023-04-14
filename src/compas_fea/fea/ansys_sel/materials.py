# Author(s): Andrew Liew (github.com/andrewliew), Marius Weber (IBK, ETHZ)

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



__all__ = [
    'Materials',
]

MPa = 10**(-6)
GPa = 10**(-9)


class Materials(object):

    def __init__(self):
        pass

    def write_materials(self):

        self.write_section('Materials')
        self.blank_line()
        
        materials = self.structure.materials
        element_set = self.structure.sets
        properties = self.structure.element_properties
        sections = self.structure.sections
        

        # Lesen der secnum
        ele_prop=[]
        
        for key_ele_prop in sorted(properties):
            ele_prop.append(key_ele_prop)

        # Lesen der benotigten Geometrien
        t_set=[]

        for key_prop in sorted(properties):
            #self.write_subsection(key_prop)

            property = properties[key_prop]
            section = sections[property.section]
            geometry = section.geometry            
            if geometry is not None:
                t_set.append(geometry.get('t', None))


        # Lesen der vorhandenen element sets
        ele_set=[]
        for key_set in sorted(self.structure.sets):
                        
            element_set = self.structure.sets[key_set]
            if element_set.type != 'node': 
                #self.write_subsection(key_set)
                ele_set.append(key_set)

        # Start Schlaufe zur Zuweisung Material zu Element sets                
        count_mat=-1
        for key, material in sorted(materials.items()):
            
            # Zuweisung Matetrial, Geometrie, Count
            E = material.E['E']
            # G = material.G
            v = material.v['v']
            p = material.p   
            t=t_set[count_mat] 
            ele_prop_name=ele_prop[count_mat]

            count_mat=count_mat+1
            G=E/(2*(1+v))
            E111=5/6*G*t    
            E221=5/6*G*t                     
            E121=0
            nn=10
            delta_h=t/nn
            
            self.write_subsection(key)
            
            mtype = material.__name__
                        
            if mtype != 'MPCStiff':
                self.write_line('!No Material properties for MPCs are needed') 
                self.blank_line()
                self.write_line('mp,ex,{0},{1}'.format(count_mat+1, E))
                self.write_line('mp,ey,{0},{1}'.format(count_mat+1, E))
                self.write_line('mp,prxy,{0},{1}'.format(count_mat+1, v))
                self.write_line('mp,dens,{0},{1}'.format(count_mat+1, p))
                self.blank_line()

                self.write_line('allsel')  
                self.write_line('cmsel,s, {0}, elem'.format(ele_set[count_mat]))  
                self.write_line('mpchg,{0},all'.format(count_mat+1))
        

            # ------------------------------------------------------------------------------------------------------
            # Ansys
            # ------------------------------------------------------------------------------------------------------

            self.blank_line()

            # Elastic
            # -------

            if mtype in ['ElasticIsotropic', 'ElasticPlastic', 'Steel', 'Concrete', 'Stiff',
                         'ConcreteSmearedCrack', 'ConcreteDamagedPlasticity']:

                self.write_line('sectype,{0} , shell'.format(count_mat+1))
                self.write_line('seccontrols,{0},{1},{2},,1,1,1'.format(E111, E221, E121))
                self.blank_line()
                self.write_line('*do,j,1,{0}'.format(nn))
                self.write_line('secdata,{0},{1},0,,,j'.format(delta_h,count_mat+1))
                self.write_line('secoffset,mid')
                self.write_line('*enddo')       

            # MPC
            # -------

            if mtype in ['MPCStiff']:

                self.write_line('! No Material properties for MPCs are needed') 


            # Compression
            # -----------

            if mtype in ['ElasticPlastic', 'Steel', 'Concrete', 'ConcreteSmearedCrack',
                         'ConcreteDamagedPlasticity']:

                raise NotImplementedError

            # Tension
            # -------

            if mtype in ['Concrete', 'ConcreteSmearedCrack']:

                raise NotImplementedError

            elif mtype == 'ConcreteDamagedPlasticity':

                raise NotImplementedError

            self.blank_line()

        self.blank_line()
        self.blank_line()

