"""An example compas_fea package use for beam elements."""

from compas_fea.cad import rhino

from compas_fea.structure import ElementProperties as Properties
from compas_fea.structure import GeneralStep
#from compas_fea.structure import LineLoad
from compas_fea.structure import PinnedDisplacement
from compas_fea.structure import PipeSection
from compas_fea.structure import PointLoad
from compas_fea.structure import RollerDisplacementXZ
from compas_fea.structure import Steel
from compas_fea.structure import Structure


__author__     = ['Andrew Liew <liew@arch.ethz.ch>']
__copyright__  = 'Copyright 2017, BLOCK Research Group - ETH Zurich'
__license__    = 'MIT License'
__email__      = 'liew@arch.ethz.ch'


# Create empty Structure object

mdl = Structure(name='beam_frame', path='C:/Temp/')

# Add nodes and elements

rhino.add_nodes_elements_from_layers(mdl, line_type='BeamElement', layers='lines')

# Add sets

rhino.add_sets_from_layers(mdl, layers=['supports', 'load_v', 'load_h', 'corners'])
rhino.add_sets_from_layers(mdl, layers='lines', explode=True)

# Add materials

mdl.add_material(Steel(name='mat_steel'))

# Add sections

mdl.add_section(PipeSection(name='sec_pipe', r=0.100, t=0.005))

# Add element properties

ep = Properties(material='mat_steel', section='sec_pipe', elsets='lines')
mdl.add_element_properties(ep, name='ep')

# Add displacements

mdl.add_displacement(PinnedDisplacement(name='disp_pins', nodes='supports'))
mdl.add_displacement(RollerDisplacementXZ(name='disp_rollers', nodes='corners'))

# Add loads

mdl.add_loads([
    PointLoad(name='load_h', nodes='load_h', x=10000),
    PointLoad(name='load_v', nodes='load_v', z=-30000)])

# Add steps

mdl.add_steps([
    GeneralStep(name='step_bc', displacements=['disp_pins', 'disp_rollers']),
    GeneralStep(name='step_loads', loads=['load_h', 'load_v'], iterations=50)])
mdl.steps_order = ['step_bc', 'step_loads']

# Structure summary

mdl.summary()

# Run and extract data

mdl.analyse_and_extract(software='opensees', fields=['u'])

# Plot

rhino.plot_data(mdl, step='step_loads', field='um', scale=50)