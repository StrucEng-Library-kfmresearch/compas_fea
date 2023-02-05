from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


# Author(s): Andrew Liew (github.com/andrewliew), Dr. Marius Weber (IBK, ETHZ)


__all__ = [
    'BCs',
]

dofs = ['x',  'y',  'z',  'xx', 'yy', 'zz']


class BCs(object):

    def __init__(self):

        pass

    def write_boundary_conditions(self):

        self.write_section('Boundary conditions')
        self.blank_line()

        sets = self.structure.sets
        steps = self.structure.steps
        displacements = self.structure.displacements

        try:

            step = steps[self.structure.steps_order[0]]

            if isinstance(step.displacements, str):
                step.displacements = [step.displacements]

            for key in step.displacements:

                nodes = displacements[key].nodes
                components = displacements[key].components
                nset = nodes if isinstance(nodes, str) else None
                selection = sets[nset].selection if isinstance(nodes, str) else nodes

                self.write_subsection(key)

 
                # ----------------------------------------------------------------------------
                # Ansys
                # ----------------------------------------------------------------------------


                self.write_line('nsel,s,,,{0}'.format(nset))
                self.blank_line()

                for c, dof in enumerate(dofs, 1):
                    if c == 1:
                            dof_ansys='ux' 
                    elif c == 2:
                            dof_ansys='uy' 
                    elif c == 3:
                            dof_ansys='uz' 
                    elif c == 4:
                            dof_ansys='rotx' 
                    elif c == 5:
                            dof_ansys='roty' 
                    elif c == 6:
                            dof_ansys='rotz'                          

                    if components[dof] is not None:
                        if nset:
                            # self.write_line('{0}, {1}, {1}, {2}'.format(nset, c, components[dof]))
                            self.write_line('d,all,{0},{1}'.format(dof_ansys,components[dof]))                                
                        else:
                            for node in sorted(selection, key=int):
                                self.write_line('{0}, {1}, {1}, {2}'.format(node + 1, c, components[dof]))
                self.blank_line()
                self.write_line('allsel')

                self.blank_line()

        except Exception:

            print('***** Error writing boundary conditions, check Step exists in structure.steps_order[0] *****')

        self.blank_line()
        self.blank_line()
