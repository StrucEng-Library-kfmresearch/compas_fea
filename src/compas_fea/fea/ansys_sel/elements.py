# Author(s): Andrew Liew (github.com/andrewliew), Marius Weber (IBK, ETHZ)

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from math import pi




__all__ = [
    'Elements',
]


ansys_data = {
    'AngleSection':       {'name': 'L',           'geometry': ['b', 'h', 't', 't']},
    'BoxSection':         {'name': 'BOX',         'geometry': ['b', 'h', 'tw', 'tf', 'tw', 'tf']},
    'CircularSection':    {'name': 'CIRC',        'geometry': ['r']},
    'ISection':           {'name': 'I',           'geometry': ['c', 'h', 'b', 'b', 'tf', 'tf', 'tw']},
    'PipeSection':        {'name': 'PIPE',        'geometry': ['r', 't']},
    'RectangularSection': {'name': 'RECTANGULAR', 'geometry': ['b', 'h']},
    'MPCSection':         {'name': 'MPC',         'geometry': None},
    'TrapezoidalSection': {'name': 'TRAPEZOID',   'geometry': ['b1', 'h', 'b2', 'c']},
    'GeneralSection':     {'name': 'GENERAL',     'geometry': ['A', 'I11', 'I12', 'I22', 'J', 'g0', 'gw']},
    'ShellSection':       {'name': None,          'geometry': ['t'], 'ORxyz': ['ORxyz'], 'XAxyz': ['XAxyz'], 'YAxyz': ['YAxyz']},
    'SolidSection':       {'name': None,          'geometry': None},
    'TrussSection':       {'name': None,          'geometry': ['A']},
}


class Elements(object):

    def __init__(self):

        pass

    def write_elements(self,structure):

        self.write_section('Elements')
        self.blank_line()
        self.write_line('allsel')
        self.blank_line()
        elements = self.structure.elements
        materials = self.structure.materials
        properties = self.structure.element_properties
        sections = self.structure.sections
        sets = self.structure.sets

        written_springs = []
        count_prop=-1 # This is only for numbering the secnum in ansys
        for key in sorted(properties):
            count_prop=count_prop+1
            self.write_subsection(key)
            self.blank_line()            
            self.write_line('secnum, {0}'.format(count_prop+1))
            self.blank_line()            
            property = properties[key]
            reinforcement = property.rebar
            elset = property.elset

            section = sections[property.section]
            stype = section.__name__
            geometry = section.geometry
            
            if geometry is not None:
                loc_coor_OR= section.loc_coords_OR
                loc_coor_XA= section.loc_coords_XA
                loc_coor_YA= section.loc_coords_YA

            material = materials.get(property.material)

            if material:
                m_index = material.index + 1

            s_index = section.index + 1

            selection = property.elements if property.elements else sets[elset].selection # Zum elemeset zugehorige Elementnummern
            
            if geometry is not None:
    
                t = geometry.get('t', None) 
                ORxyz = loc_coor_OR.get('ORxyz', None)
                XAxyz = loc_coor_XA.get('XAxyz', None)
                YAxyz = loc_coor_YA.get('YAxyz', None)                    
                A = geometry.get('A', None)
                J = geometry.get('J', None)
                Ixx = geometry.get('Ixx', None)
                Iyy = geometry.get('Iyy', None)
                E = material.E.get('E', None)
                G = material.G.get('G', None)
                v = material.v['v']
                p = material.p   

            for select in selection:

                element = elements[select]
                nodes = [str(i + 1) for i in element.nodes]
                no = len(nodes)
                n = select + 1
                
                #ex = element.axes.get('ex', None)
                #ey = element.axes.get('ey', None)
                #ez = element.axes.get('ez', None)

                # =====================================================================================================
                # =====================================================================================================
                # SOLID
                # =====================================================================================================
                # =====================================================================================================

                if stype == 'SolidSection':
                    raise NotImplementedError

                # =====================================================================================================
                # =====================================================================================================
                # SHELL
                # =====================================================================================================
                # =====================================================================================================

                elif stype == 'ShellSection':


                        
                    e = 'element_{0}'.format(select)
                    self.write_line('! ----- Start Element {0} ----- '.format(n))
                    
                    # Create element
                    self.write_line('et, {0}, shell181'.format(n))
                    self.write_line('type, {0}'.format(n))     
                    self.write_line('en,{0},{1}'.format(n, ','.join(nodes)))
                    self.write_line('keyopt, {0}, 1,0'.format(n))
                    self.write_line('keyopt, {0}, 3,2'.format(n))
                    self.write_line('keyopt, {0}, 8,2'.format(n))

                    # loc coor system for the element
                    self.write_line('esel,s,elem,,{0}'.format(n))
                    self.write_line('k,1,{0},{1},{2}'.format(ORxyz[0],ORxyz[1],ORxyz[2]))   #set kp origin          
                    self.write_line('k,2,{0},{1},{2}'.format(XAxyz[0],XAxyz[1],XAxyz[2]))   #set kp in x direction
                    self.write_line('k,3,{0},{1},{2}'.format(YAxyz[0],YAxyz[1],YAxyz[2]))   #set kp in y direction                    
                    No_loc_system=10+n
                    self.write_line('cskp,{0},0,1,2,3'.format(No_loc_system)) # Set and actiave local coordiante system                   
                    self.write_line('emodif,all,esys,{}'.format(No_loc_system))
                    self.write_line('kdele,1,3,1') 
                    self.write_line('allsel') 
                    self.write_line('csys,0') # set to the orignal coordiante system                   
    


                    if reinforcement:
                        raise NotImplementedError


                # =====================================================================================================
                # =====================================================================================================
                # TRUSS
                # =====================================================================================================
                # =====================================================================================================
                
                elif stype == 'TrussSection':
                    raise NotImplementedError

                # =====================================================================================================
                # =====================================================================================================
                # SPRING
                # =====================================================================================================
                # =====================================================================================================

                elif stype == 'SpringSection':

                    raise NotImplementedError

                # =====================================================================================================
                # =====================================================================================================
                # MASS
                # =====================================================================================================
                # =====================================================================================================

                elif stype == 'MassSection':

                    raise NotImplementedError

                # =====================================================================================================
                # =====================================================================================================
                # MPC 
                # =====================================================================================================
                # =====================================================================================================

                elif stype == 'MPCSection':

                    self.write_line('! ----- Start Element {0} ----- '.format(n))                                      
                    
                    # Create element
                    self.write_line('et, {0}, mpc184'.format(n))                    
                    self.write_line('type, {0}'.format(n))             
                    self.write_line('en,{0},{1}'.format(n, ','.join(nodes)))
                    self.write_line('keyopt, {0}, 1,1'.format(n))

                # =====================================================================================================
                # =====================================================================================================
                # BEAM
                # =====================================================================================================
                # =====================================================================================================

                else:
                                        
                    raise NotImplementedError


            self.blank_line()
            self.blank_line()
