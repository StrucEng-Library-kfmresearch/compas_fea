from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


# Author(s): Andrew Liew (github.com/andrewliew)


__all__ = [
    'Sets',
]


class Sets(object):

    def __init__(self):

        pass

    def write_node_sets(self):

        self.write_section('Node sets')
        self.blank_line()

        for key in sorted(self.structure.sets):

            node_set = self.structure.sets[key]

            if node_set.type == 'node':
                self.write_node_set(key, node_set)
                self.blank_line()
                self.blank_line()

    def write_node_set(self, key, node_set):

        self.write_subsection(key)

        header = {'ansys_sel':   'allsel'}

        self.write_line(header[self.software])
        self.blank_line()

        nodes = [i + 1 for i in node_set.selection]
        
        for i in range(0, len(nodes), 1):
            
            for j in nodes[i:i + 1]:
                node_string=str(j) 
                self.write_line('nsel,u,node,, {0}'.format(node_string))
        self.blank_line()
        self.write_line('nsel,inve')  
        self.write_line('CM, {0},node'.format(key))


    def write_element_sets(self):

        self.write_section('Element sets')
        self.blank_line()

        for key in sorted(self.structure.sets):

            element_set = self.structure.sets[key]

            if element_set.type != 'node':
                self.write_element_set(key, element_set)
                self.blank_line()
                self.blank_line()

    def write_element_set(self, key, element_set):

        stype = element_set.type

        self.write_subsection(key)

        if stype in ['element', 'surface_node']:

            if stype == 'element':
                self.write_line('allsel') #('*ELSET, ELSET={0}'.format(key))

            elif stype == 'surface_node':
                raise NotImplementedError

            self.blank_line()

            selection = [i + 1 for i in sorted(element_set.selection)]

            # Set for elements
            for i in range(0, len(selection), 1):

                for j in selection[i:i + 1]:
                    element_string=str(j) 
                    self.write_line('esel,u,elem,, {0}'.format(element_string))
            self.blank_line()
            self.write_line('esel,inve')  
            self.write_line('CM, {0},elem'.format(key))     


        if stype == 'surface_element':

            raise NotImplementedError
