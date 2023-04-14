# Author(s): Andrew Liew (github.com/andrewliew), Marius Weber (IBK, ETHZ)

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function





__all__ = [
    'Heading',
]


class Heading(object):

    def __init__(self):
        pass

    def write_heading(self):
        header = {'ansys_sel':   'finish\n/clear\n/prep7\n' }
        
        self.write_section('Heading')
        self.blank_line()
        self.write_line(header[self.software])
        self.blank_line()
        self.blank_line()
        self.blank_line()

