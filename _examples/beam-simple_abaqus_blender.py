"""An example compas_fea package use for beam elements."""

from math import pi

from compas_fea.fea.abaq import abaq
from compas_fea.cad import blender

from compas_fea.structure import CircularSection
from compas_fea.structure import ElasticIsotropic
from compas_fea.structure import ElementProperties
from compas_fea.structure import GeneralDisplacement
from compas_fea.structure import GeneralStep
from compas_fea.structure import PointLoad
from compas_fea.structure import Structure

from compas.cad.blender.geometry import add_empty
from compas.cad.blender.helpers import network_from_bmesh
from compas.cad.blender.utilities import draw_bmesh
from compas.cad.blender.utilities import clear_layers


__author__     = ['Andrew Liew <liew@arch.ethz.ch>']
__copyright__  = 'Copyright 2017, BLOCK Research Group - ETH Zurich'
__license__    = 'MIT License'
__email__      = 'liew@arch.ethz.ch'


# Folders and Structure name

name = 'beam-simple'
path = '/home/al/Temp/'

# Create an empty Structure object named mdl

mdl = Structure()

# Clear layers

clear_layers([0, 1, 2, 3])

# Create beam

L = 1.0
nx = 100
x = [i * L / nx for i in range(nx + 1)]
vertices = [[xi, 0, 0] for xi in x]
edges = [[i, i + 1] for i in range(nx)]
bmesh = draw_bmesh(name='beam', vertices=vertices, edges=edges, layer=0)
network = network_from_bmesh(bmesh)
add_empty(location=[0, 0, 0], layer=1, radius=0.05)
add_empty(location=[L, 0, 0], layer=2, radius=0.05)

# Weights

spacing = 5
weight = -1.0
for xi in x[spacing::spacing]:
    add_empty(type='SPHERE', location=[xi, 0, 0], layer=3, radius=0.005)

# Add beam element and geometry to Structure

mdl.add_nodes_elements_from_network(network=network, element_type='BeamElement')

# Add node and element sets to the Structure object

blender.add_elset_from_bmeshes(mdl, layer=0, name='elset_lines', explode=True)
blender.add_nset_from_objects(mdl, layer=1, name='nset_left')
blender.add_nset_from_objects(mdl, layer=2, name='nset_right')
blender.add_nset_from_objects(mdl, layer=3, name='nset_weights')

# Add materials

E = 20 * 10**9
mdl.add_material(ElasticIsotropic(name='mat_elastic', E=E, v=0.3, p=1500))

# Add sections

nodes, elements, arclengths, L = blender.ordered_lines(mdl, network=network, layer=1)
for c, i in enumerate(arclengths):
    ri = (((i / L) - 0.5)**2 + 0.4) * 0.01
    sname = 'sec_rod_{0}'.format(elements[c])
    mdl.add_section(CircularSection(name=sname, r=ri))
    ep = ElementProperties(material='mat_elastic', section=sname, elsets='element_{0}'.format(elements[c]))
    mdl.add_element_properties(ep, name='ep_{0}'.format(elements[c]))

# Add loads

mdl.add_load(PointLoad(name='load_weights', nodes='nset_weights', z=weight))

# Add displacements

deg = pi / 180
mdl.add_displacement(GeneralDisplacement(name='disp_left', nodes='nset_left', x=0, y=0, z=0, xx=0, zz=0))
mdl.add_displacement(GeneralDisplacement(name='disp_right', nodes='nset_right', y=0, z=-0.2, xx=0, zz=0))

# Add steps

mdl.add_step(GeneralStep(name='step_bc', displacements=['disp_left', 'disp_right']))
mdl.add_step(GeneralStep(name='step_load', loads=['load_weights']))
mdl.set_steps_order(['step_bc', 'step_load'])

# Structure summary

mdl.summary()

# Generate .inp file

fnm = '{0}{1}.inp'.format(path, name)
abaq.inp_generate(mdl, filename=fnm)

# Run and extract data

exe = '/home/al/abaqus/Commands/abaqus cae '
mdl.analyse(path=path, name=name, software='abaqus', exe=exe, fields='U,UR,SF,SM')

# Plot displacements

blender.plot_data(mdl, path, name, step='step_load', field='U', component='magnitude', radius=0.01, layer=4)
blender.plot_data(mdl, path, name, step='step_load', field='UR', component='UR2', radius=0.01, layer=5)

# Plot section forces/moments

blender.plot_data(mdl, path, name, step='step_load', field='SF', component='SF1', radius=0.01, layer=6)
blender.plot_data(mdl, path, name, step='step_load', field='SF', component='SF3', radius=0.01, layer=7)
blender.plot_data(mdl, path, name, step='step_load', field='SM', component='SM2', radius=0.01, layer=8)
