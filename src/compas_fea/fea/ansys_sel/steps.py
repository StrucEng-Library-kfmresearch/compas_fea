# Author(s): Andrew Liew (github.com/andrewliew), Marius Weber (IBK, ETHZ)

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import re

import os
import json


__all__ = [
    'Steps',
]

dofs = ['x', 'y', 'z', 'xx', 'yy', 'zz']


class Steps(object):

    def __init__(self):
        pass

    def write_steps(self):

        self.write_section('Steps')
        self.blank_line()

        displacements = self.structure.displacements
        loads = self.structure.loads
        steps = self.structure.steps
        sets = self.structure.sets
        fields = self.fields

        # temp folder

        temp = '{0}{1}/'.format(self.structure.path, self.structure.name)

        try:
            os.stat(temp)
            for file in os.listdir(temp):
                os.remove(os.path.join(temp, file))
        except Exception:
            os.mkdir(temp)

        # Steps

        for key in self.structure.steps_order[1:]:

            step = steps[key]
            stype = step.__name__
            s_index = step.index
            factor = getattr(step, 'factor', 1)
            increments = getattr(step, 'increments', 100)
            iterations = getattr(step, 'iterations', 100)
            tolerance = getattr(step, 'tolerance', None)
            method = getattr(step, 'type')
            modes = getattr(step, 'modes', None)
            modify = getattr(step, 'modify', None)
            nlgeom = 'YES' if getattr(step, 'nlgeom', None) else 'NO'
            op = 'MOD' if modify else 'NEW'

            # =====================================================================================================
            # =====================================================================================================
            # HEADER
            # =====================================================================================================
            # =====================================================================================================
            self.write_subsection(key)
                


            # =====================================================================================================
            # =====================================================================================================
            # LOADS
            # =====================================================================================================
            # =====================================================================================================

            if getattr(step, 'loads', None):

                if isinstance(step.loads, str):
                    step.loads = [step.loads]

                for k in step.loads:

                    self.write_subsection(k)

                    load = loads[k]
                    ltype = load.__name__
                    com = getattr(load, 'components', None)
                    axes = getattr(load, 'axes', None)
                    nodes = getattr(load, 'nodes', None)
                    fact = factor.get(k, 1.0) if isinstance(factor, dict) else factor

                    if com:
                        gx = com.get('x', 0)
                        gy = com.get('y', 0)
                        gz = com.get('z', 0)

                    if isinstance(nodes, str):
                        nodes = [nodes]

                    if isinstance(load.elements, str):
                        elements = [load.elements]
                    else:
                        elements = load.elements

                    # -------------------------------------------------------------------------------------------------
                    # Ansys
                    # -------------------------------------------------------------------------------------------------

                    # PointLoad
                    # ---------

                    if ltype == 'PointLoad':



                        #self.write_line('*CLOAD, OP={0}'.format(op))
                        #self.blank_line()

                        for node in nodes:

                            ni = node if isinstance(node, str) else node + 1
                    
                            self.write_line('nsel,s,,,{0}'.format(ni))                            
                            self.blank_line()

                            for c, dof in enumerate(dofs, 1):
                                if c == 1:
                                    dof_ansys='fx' 
                                elif c == 2:
                                    dof_ansys='fy' 
                                elif c == 3:
                                    dof_ansys='fz' 
                                elif c == 4:
                                    dof_ansys='mx' 
                                elif c == 5:
                                    dof_ansys='my' 
                                elif c == 6:
                                    dof_ansys='mz'   

                                if com[dof]:                                                                                                                            
                                    self.write_line('F,all,{0},{1}'.format(dof_ansys, com[dof] * fact))
                            
                            self.write_line('allsel')

                        # AreaLoad
                        # --------

                    elif ltype == 'AreaLoad':

                        for k in elements:

                            if com['z']:
                                
                                self.write_line('allsel')  
                                #self.write_line('cmsel,s, {0}, elem'.format(k))  
                                self.write_line('esel,s,elem,,{0}'.format(k))  
                                self.write_line('sfcum,pres,add,,,')
                                self.write_line('sfe,all,1,pres,1,{0}'.format(fact * com['z']))
                                self.write_line('allsel')  
                                self.blank_line()

                            elif com['x']:
                                
                                self.write_line('allsel')  
                                #self.write_line('cmsel,s, {0}, elem'.format(k))  
                                self.write_line('esel,s,elem,,{0}'.format(k))  
                                self.write_line('sfcum,pres,add,,,')
                                self.write_line('sfe,all,6,pres,1,{0}'.format(fact * com['x']))
                                self.write_line('allsel')  
                                self.blank_line()
                            
                            elif com['y']:
                                
                                self.write_line('allsel')  
                                #self.write_line('cmsel,s, {0}, elem'.format(k))  
                                self.write_line('esel,s,elem,,{0}'.format(k))  
                                self.write_line('sfcum,pres,add,,,')
                                self.write_line('sfe,all,3,pres,1,{0}'.format(fact * com['y']))
                                self.write_line('allsel')  
                                self.blank_line()                                
                            else:
                                raise NotImplementedError                                                                


                    # PointLoads
                    # ----------

                    elif ltype == 'PointLoads':
                        self.write_line('*CLOAD, OP={0}'.format(op))
                        self.blank_line()

                        for node, coms in com.items():
                            for ci, value in coms.items():
                                index = dofs.index(ci) + 1
                                self.write_line('{0}, {1}, {2}'.format(node + 1, index, value * fact))

                    # Gravity
                    # -------

                    elif ltype == 'GravityLoad':

                        for k in elements:

                            #self.write_line('*DLOAD, OP={0}'.format(op))
                            self.blank_line()
                            self.write_line('acel,{0}, {1}, {2}'.format(9.81*gx, 9.81*gy, 9.81*gz))
                            #self.write_line('{0}, GRAV, {1}, {2}, {3}, {4}'.format(k, -9.81 * fact, gx, gy, gz))
                            self.blank_line()

                    # TributaryLoad
                    # -------------

                    elif ltype == 'TributaryLoad':

                        self.write_line('*CLOAD, OP={0}'.format(op))
                        self.blank_line()

                        for node in sorted(com, key=int):

                            ni = node + 1

                            for ci, dof in enumerate(dofs[:3], 1):
                                if com[node][dof]:
                                    self.write_line('{0}, {1}, {2}'.format(ni, ci, com[node][dof] * fact))

                    # LineLoad
                    # --------

                    elif ltype == 'LineLoad':

                        for k in elements:

                            self.write_line('*DLOAD, OP={0}'.format(op))
                            self.blank_line()

                            if axes == 'global':

                                for dof in dofs[:3]:
                                    if com[dof]:
                                        self.write_line('{0}, P{1}, {2}'.format(k, dof.upper(), fact * com[dof]))

                            elif axes == 'local':

                                if com['x']:
                                    self.write_line('{0}, P1, {1}'.format(k, fact * com['x']))

                                if com['y']:
                                    self.write_line('{0}, P2, {1}'.format(k, fact * com['y']))

                    # Prestress
                    # ---------

                    elif ltype == 'PrestressLoad':

                        for k in elements:

                            stresses = ''

                            if com['sxx']:
                                stresses += str(com['sxx'] * fact)

                            self.write_line('*INITIAL CONDITIONS, TYPE=STRESS')
                            self.blank_line()
                            self.write_line('{0}, {1}'.format(k, stresses))

            # -------------------------------------------------------------------------------------------------
            # Ansys Solver Block
            # -------------------------------------------------------------------------------------------------

            if stype in ['GeneralStep', 'BucklingStep', 'ModalStep']:

                if stype == 'ModalStep':

                    raise NotImplementedError

                else:

                    p = ', PERTURBATION' if stype == 'BucklingStep' else ''
                    self.write_line('/solu')
                    self.write_line('cnvtol,F,,0.9')
                    self.write_line('cnvtol,U,,0.9')
                    self.write_line('cnvtol,M,-1,3')
                    self.write_line('autots,1')
                    self.write_line('nsubst,{0},{0},{0}'.format(increments))
                    ansys_step=re.sub('step_','',key)
                    ansys_step_string=int(ansys_step)-1
                    self.write_line('time, {0}'.format(ansys_step_string))
                    self.write_line('Nropt,Full,,on')
                    self.write_line('NLGEOM,off')

                    if stype == 'BucklingStep':
                        raise NotImplementedError

                self.blank_line()

                self.blank_line()
                self.blank_line()
                self.blank_line()

            # =====================================================================================================
            # =====================================================================================================
            # DISPLACEMENTS
            # =====================================================================================================
            # =====================================================================================================

            if getattr(step, 'displacements', None):

                if isinstance(step.displacements, str):
                    step.displacements = [step.displacements]

                for k in step.displacements:

                    displacement = displacements[k]
                    com = displacement.components
                    nodes = displacement.nodes

                    if isinstance(nodes, str):
                        nodes = [nodes]

                    fact = factor.get(k, 1.0) if isinstance(factor, dict) else factor

                    # -------------------------------------------------------------------------------------------------
                    # Ansys
                    # -------------------------------------------------------------------------------------------------

                    if stype not in ['ModalStep', 'BucklingStep']:

                        self.write_line('*BOUNDARY')
                        self.blank_line()

                        for node in nodes:

                            ni = node if isinstance(node, str) else node + 1

                            for c, dof in enumerate(dofs, 1):
                                if com[dof] is not None:
                                    self.write_line('{0}, {1}, {1}, {2}'.format(ni, c, com[dof] * fact))

                        self.blank_line()
                        self.blank_line()


            # =====================================================================================================
            # =====================================================================================================
            # OUTPUT
            # =====================================================================================================
            # =====================================================================================================

            self.write_subsection('Output')

            # -------------------------------------------------------------------------------------------------
            # Ansys
            # -------------------------------------------------------------------------------------------------

           
            node_fields = ['rf', 'rm', 'u', 'ur', 'cf', 'cm']
            element_fields = ['sf', 'sm', 'sk', 'se', 's', 'e', 'pe', 'rbfor', 'ctf']

            if 'spf' in fields:
                fields[fields.index('spf')] = 'ctf'

            self.write_line('outres,erase')
            self.write_line('outres,nsol,last')
            self.write_line('outres,rsol,last')
            self.write_line('outres,esol,last')
            self.write_line('outres,svar,all')

            if 'rbfor' in fields:
                self.write_line('*ELEMENT OUTPUT, REBAR')
                self.write_line('RBFOR')

            self.blank_line()
            self.write_line('solve')
            self.blank_line()
            self.blank_line()
