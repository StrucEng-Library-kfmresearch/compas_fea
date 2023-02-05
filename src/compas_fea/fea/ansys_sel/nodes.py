from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


# Author(s): Andrew Liew (github.com/andrewliew), Marius Weber (IBK, ETHZ)


__all__ = [
    'Nodes',
]


class Nodes(object):

    def __init__(self):
        pass

    def write_nodes(self):

        header = {'ansys_sel':   '!'}

        self.prefix = {'ansys_sel':   ''}

        self.write_section('Nodes')
        self.write_line(header[self.software])

        for key in sorted(self.structure.nodes, key=int):

            self.write_node(key)

        self.blank_line()
        self.blank_line()

    def write_node(self, key):

        prefix = self.prefix[self.software]
        spacer = self.spacer[self.software]
        x, y, z = self.structure.node_xyz(key)

        line_1 = '{0}{1}{2}{3:.3f}{2}{4:.3f}{2}{5:.3f}'.format('n,', key + 1, spacer, x, y, z)
        self.write_line(line_1)


    def write_mass(self, key):

        raise NotImplementedError


