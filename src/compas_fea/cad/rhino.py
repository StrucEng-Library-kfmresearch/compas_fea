# Author(s): Compas/Compas FEA Team, Marius  Weber (ETHZ, HSLU T&A)

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import json

import compas
if compas.RHINO:
    from compas_rhino.geometry import RhinoMesh

from compas.datastructures.mesh import Mesh
from compas.datastructures import Network
from compas.geometry import Frame, Transformation, Vector, matrix_from_axis_and_angle
from compas.geometry import add_vectors
from compas.geometry import cross_vectors
from compas.geometry import normalize_vector
from compas.geometry import length_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas_fea.structure import Structure
from time import time
import math

from compas_fea.utilities import colorbar
from compas_fea.utilities import extrude_mesh
from compas_fea.utilities import network_order
from compas.rpc import Proxy
import statistics

if not compas.IPY:
    from compas_fea.utilities import meshing
    from compas_fea.utilities import functions

else:
    from compas.rpc import Proxy
    functions = Proxy('compas_fea.utilities.functions')
    meshing = Proxy('compas_fea.utilities.meshing')
    

if compas.RHINO:
    import rhinoscriptsyntax as rs



__all__ = [
    'add_element_set',
    'add_node_set',
    'add_nodes_elements_from_layers',
    'add_sets_from_layers',
    'add_tets_from_mesh',
    'discretise_mesh',
    'mesh_extrude',
    'network_from_lines',
    'ordered_network',
    'plot_reaction_forces',
    'plot_concentrated_forces',
    'plot_mode_shapes',
    'plot_volmesh',
    'plot_axes',
    'plot_data',
    'plot_principal_stresses',
    'plot_principal_strains',
    'plot_steel_stresses',
    'plot_steel_shear',
    'plot_voxels',
    'weld_meshes_from_layer',
]


def add_element_set(structure, guids, name):
    """
    Adds element set information from Rhino curve and mesh guids.

    Parameters
    ----------
    structure : obj
        Structure object to update.
    guids : list
        Rhino curve and Rhino mesh guids.
    name : str
        Name of the new element set.

    Returns
    -------
    None

    Notes
    -----
    - Meshes representing solids must have 'solid' in their name.
    """

    elements = []

    for guid in guids:

        if rs.IsCurve(guid):

            sp = structure.check_node_exists(rs.CurveStartPoint(guid))
            ep = structure.check_node_exists(rs.CurveEndPoint(guid))
            element = structure.check_element_exists([sp, ep])
            if element is not None:
                elements.append(element)

        elif rs.IsMesh(guid):

            vertices = rs.MeshVertices(guid)
            faces = rs.MeshFaceVertices(guid)

            if rs.ObjectName(guid) and ('solid' in rs.ObjectName(guid)):
                nodes = [structure.check_node_exists(i) for i in vertices]
                element = structure.check_element_exists(nodes)
                if element is not None:
                    elements.append(element)
            else:
                for face in faces:
                    nodes = [structure.check_node_exists(vertices[i]) for i in face]
                    if nodes[2] == nodes[3]:
                        del nodes[-1]
                    element = structure.check_element_exists(nodes)
                    if element is not None:
                        elements.append(element)

    structure.add_set(name=name, type='element', selection=elements)


def add_node_set(structure, guids, name):
    """
    Adds node set information from Rhino point guids.

    Parameters
    ----------
    structure : obj
        Structure object to update.
    guids : list
        Rhino point guids.
    name : str
        Name of the new node set.

    Returns
    -------
    None

    """

    nodes = []

    for guid in guids:

        if rs.IsPoint(guid):

            node = structure.check_node_exists(rs.PointCoordinates(guid))
            if node is not None:
                nodes.append(node)

    structure.add_set(name=name, type='node', selection=nodes)


def add_nodes_elements_from_layers(structure, layers, line_type=None, mesh_type=None, thermal=False, pA=None, pL=None):
    """
    Adds node and element data from Rhino layers to the Structure object.

    Parameters
    ----------
    structure : obj
        Structure object to update.
    layers : list
        Layer string names to extract nodes and elements.
    line_type : str
        Element type for line objects.
    mesh_type : str
        Element type for mesh objects.
    thermal : bool
        Thermal properties on or off.
    pA : float
        Mass area density [kg/m2].
    pL : float
        Mass length density [kg/m].

    Returns
    -------
    list
        Node keys that were added to the Structure.
    list
        Element keys that were added to the Structure.

    """

    if isinstance(layers, str):
        layers = [layers]

    added_nodes = set()
    added_elements = set()

    for layer in layers:

        elset = set()

        for guid in rs.ObjectsByLayer(layer):

            if line_type and rs.IsCurve(guid):

                sp_xyz = rs.CurveStartPoint(guid)
                ep_xyz = rs.CurveEndPoint(guid)
                ez = subtract_vectors(ep_xyz, sp_xyz)
                L = length_vector(ez)
                m = 0.5 * L * pL if pL else None

                sp = structure.add_node(xyz=sp_xyz, mass=m)
                ep = structure.add_node(xyz=ep_xyz, mass=m)
                added_nodes.add(sp)
                added_nodes.add(ep)

                try:
                    name = rs.ObjectName(guid).replace("'", '"')

                    if name[0] in ['_', '^']:
                        name = name[1:]

                    dic = json.loads(name)
                    ex = dic.get('ex', None)
                    ey = dic.get('ey', None)

                    if ex and not ey:
                        ey = cross_vectors(ex, ez)

                except Exception:
                    ex = None
                    ey = None

                axes = {'ex': ex, 'ey': ey, 'ez': ez}

                ekey = structure.add_element(nodes=[sp, ep], type=line_type, thermal=thermal, axes=axes)

                if (line_type == 'BeamElement') and (ex is None):

                    if (ez[0] == 0) and (ez[1] == 0):

                        print('***** WARNING: vertical BeamElement with no ex axis, element {0} *****'.format(ekey))

                if ekey is not None:
                    added_elements.add(ekey)
                    elset.add(ekey)

            elif mesh_type and rs.IsMesh(guid):

                mesh = RhinoMesh.from_guid(guid).to_compas()

                vertices = rs.MeshVertices(guid)
                nodes = []
                masses = []

                for c, vertex in enumerate(vertices):
                    m = mesh.vertex_area(c) * pA if pA else None
                    masses.append(m)
                    nodes.append(structure.add_node(xyz=vertex, mass=m))

                added_nodes.update(nodes)

                if mesh_type in ['HexahedronElement', 'TetrahedronElement', 'SolidElement', 'PentahedronElement']:
                    ekey = structure.add_element(nodes=nodes, type=mesh_type, thermal=thermal)

                    if ekey is not None:
                        added_elements.add(ekey)
                        elset.add(ekey)

                elif mesh_type == 'MassElement':
                    node_iterator = 0
                    for node in nodes:
                        # structure.nodes[node].mass
                        ekey = structure.add_element(nodes=[node], type=mesh_type,
                                                     thermal=thermal, mass=masses[node_iterator])
                        node_iterator += 1
                        if ekey is not None:
                            added_elements.add(ekey)
                            elset.add(ekey)

                else:

                    try:
                        name = rs.ObjectName(guid).replace("'", '"')

                        if name[0] in ['_', '^']:
                            name = name[1:]

                        dic = json.loads(name)
                        ex = dic.get('ex', None)
                        ey = dic.get('ey', None)
                        ez = dic.get('ez', None)

                        if (ex and ey) and (not ez):
                            ez = cross_vectors(ex, ey)

                    except Exception:
                        ex = None
                        ey = None
                        ez = None

                    axes = {'ex': ex, 'ey': ey, 'ez': ez}

                    for face in rs.MeshFaceVertices(guid):

                        nodes = [structure.check_node_exists(vertices[i]) for i in face]
                        if nodes[-1] == nodes[-2]:
                            del nodes[-1]

                        ekey = structure.add_element(nodes=nodes, type=mesh_type, thermal=thermal, axes=axes)
                        if ekey is not None:
                            added_elements.add(ekey)
                            elset.add(ekey)

        structure.add_set(name=layer, type='element', selection=list(elset))

    return list(added_nodes), list(added_elements)


def add_sets_from_layers(structure, layers):
    """
    Add node and element sets to the Structure object from Rhino layers.

    Parameters
    ----------
    structure : obj
        Structure object to update.
    layers : list
        List of layer string names to take objects from.

    Returns
    -------
    None

    Notes
    -----
    - Layers should exclusively contain nodes or elements.
    - Mixed elements, e.g. lines and meshes, are allowed on a layer.
    - Sets will inherit the layer names as their set name.

    """

    if isinstance(layers, str):
        layers = [layers]

    for layer in layers:
        guids = rs.ObjectsByLayer(layer)

        if guids:
            name = layer.split('::')[-1] if '::' in layer else layer
            check_points = [rs.IsPoint(guid) for guid in guids]

            if all(check_points):
                add_node_set(structure=structure, guids=guids, name=name)
            elif not any(check_points):
                add_element_set(structure=structure, guids=guids, name=name)
            else:
                print('***** Layer {0} contained a mixture of points and elements, set not created *****'.format(name))


def add_tets_from_mesh(structure, name, mesh, draw_tets=False, volume=None, thermal=False):
    """
    Adds tetrahedron elements from a mesh in Rhino to the Structure object.

    Parameters
    ----------
    structure : obj
        Structure object to update.
    name : str
        Name for the element set of tetrahedrons.
    mesh : guid
        The mesh in Rhino representing the outer surface.
    draw_tets : str
        Layer to draw tetrahedrons on.
    volume : float
        Maximum volume for each tet.
    thermal : bool
        Thermal properties on or off.

    Returns
    -------
    None

    """

    rhinomesh = RhinoMesh.from_guid(mesh)
    vertices = rhinomesh.vertices
    faces = [face[:3] for face in rhinomesh.faces]

    try:
        tets_points, tets_elements = meshing.tets_from_vertices_faces(vertices=vertices, faces=faces, volume=volume)

        for point in tets_points:
            structure.add_node(point)

        ekeys = []

        for element in tets_elements:

            nodes = [structure.check_node_exists(tets_points[i]) for i in element]
            ekey = structure.add_element(nodes=nodes, type='TetrahedronElement', thermal=thermal)
            ekeys.append(ekey)

        structure.add_set(name=name, type='element', selection=ekeys)

        if draw_tets:

            rs.EnableRedraw(False)
            rs.DeleteObjects(rs.ObjectsByLayer(draw_tets))
            rs.CurrentLayer(draw_tets)

            tet_faces = [[0, 2, 1, 1], [1, 2, 3, 3], [1, 3, 0, 0], [0, 3, 2, 2]]

            for i, points in enumerate(tets_elements):

                xyz = [tets_points[j] for j in points]
                rs.AddMesh(vertices=xyz, face_vertices=tet_faces)

            rs.EnableRedraw(True)

        print('***** MeshPy (TetGen) successfull *****')

    except Exception:

        print('***** Error using MeshPy (TetGen) or drawing Tets *****')


def discretise_mesh(mesh, layer, target, min_angle=15, factor=1):
    """
    Discretise a mesh from an input triangulated coarse mesh into small denser meshes.

    Parameters
    ----------
    mesh : guid
        The guid of the Rhino input mesh.
    layer : str
        Layer name to draw results.
    target : float
        Target length of each triangle.
    min_angle : float
        Minimum internal angle of triangles.
    factor : float
        Factor on the maximum area of each triangle.

    Returns
    -------
    None

    """

    rhinomesh = RhinoMesh.from_guid(mesh)
    vertices = rhinomesh.vertices
    faces = [face[:3] for face in rhinomesh.faces]

    try:

        points, tris = meshing.discretise_faces(vertices=vertices, faces=faces,
                                                target=target, min_angle=min_angle, factor=factor)

        rs.CurrentLayer(rs.AddLayer(layer))
        rs.DeleteObjects(rs.ObjectsByLayer(layer))
        rs.EnableRedraw(False)

        for pts, tri in zip(points, tris):
            mesh_faces = []

            for i in tri:
                face_ = i + [i[-1]]
                mesh_faces.append(face_)
            rs.AddMesh(pts, mesh_faces)

        rs.EnableRedraw(True)

    except Exception:

        print('***** Error using MeshPy (Triangle) or drawing faces *****')


def mesh_extrude(structure, guid, layers, thickness, mesh_name='', links_name='', blocks_name='', points_name='',
                 plot_mesh=False, plot_links=False, plot_blocks=False, plot_points=False):
    """
    Extrudes a Rhino mesh and adds/creates elements.

    Parameters
    ----------
    structure : obj
        Structure object to update.
    guid : guid
        Rhino mesh guid.
    layers : int
        Number of layers.
    thickness : float
        Layer thickness.
    mesh_name : str
        Name of set for mesh on final surface.
    links_name : str
        Name of set for adding links along extrusion.
    blocks_name : str
        Name of set for solid elements.
    points_name : str
        Name of aded points.
    plot_mesh : bool
        Plot outer mesh.
    plot_links : bool
        Plot links.
    plot_blocks : bool
        Plot blocks.
    plot_points : bool
        Plot end points.

    Returns
    -------
    None

    Notes
    -----
    - Extrusion is along the mesh vertex normals.

    """

    mesh = RhinoMesh.from_guid(guid).to_compas(cls=Mesh)
    extrude_mesh(structure=structure, mesh=mesh, layers=layers, thickness=thickness, mesh_name=mesh_name,
                 links_name=links_name, blocks_name=blocks_name)

    block_faces = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]
    xyz = structure.nodes_xyz()

    rs.EnableRedraw(False)

    if plot_blocks:

        rs.CurrentLayer(rs.AddLayer(blocks_name))
        rs.DeleteObjects(rs.ObjectsByLayer(blocks_name))

        for i in structure.sets[blocks_name]['selection']:
            nodes = structure.elements[i].nodes
            xyz = structure.nodes_xyz(nodes)
            rs.AddMesh(xyz, block_faces)

    if plot_mesh:

        rs.CurrentLayer(rs.AddLayer(mesh_name))
        rs.DeleteObjects(rs.ObjectsByLayer(mesh_name))

        faces = []
        for i in structure.sets[mesh_name]['selection']:
            enodes = structure.elements[i].nodes
            if len(enodes) == 3:
                enodes.append(enodes[-1])
            faces.append(enodes)
        rs.AddMesh(xyz, faces)

    if plot_links:

        rs.CurrentLayer(rs.AddLayer(links_name))
        rs.DeleteObjects(rs.ObjectsByLayer(links_name))

        if plot_points:
            rs.CurrentLayer(rs.AddLayer(points_name))
            rs.DeleteObjects(rs.ObjectsByLayer(points_name))

        for i in structure.sets[links_name]['selection']:

            nodes = structure.elements[i].nodes
            xyz = structure.nodes_xyz(nodes)
            rs.CurrentLayer(links_name)
            rs.AddLine(xyz[0], xyz[1])

            if plot_points:
                rs.CurrentLayer(points_name)
                rs.AddPoint(xyz[1])

    rs.EnableRedraw(True)
    rs.CurrentLayer(rs.AddLayer('Default'))


def network_from_lines(guids=[], layer=None):
    """
    Creates a Network datastructure object from a list of Rhino curve guids.

    Parameters
    ----------
    guids : list
        guids of the Rhino curves to be made into a Network.
    layer : str
        Layer to grab line guids from.

    Returns
    -------
    obj
        Network datastructure object.

    """

    if layer:
        guids = rs.ObjectsByLayer(layer)
    lines = [[rs.CurveStartPoint(guid), rs.CurveEndPoint(guid)] for guid in guids if rs.IsCurve(guid)]

    return Network.from_lines(lines)


def ordered_network(structure, network, layer):
    """
    Extract vertex and edge orders from a Network for a given start-point.

    Parameters
    ----------
    structure : obj
        Structure object.
    network : obj
        Network Datastructure object.
    layer : str
        Layer to extract start-point (Rhino point).
    Returns
    -------
    list
        Ordered nodes for the Structure.
    list
        Ordered elements for the Structure.
    list
        Cumulative length at element mid-points.
    float
        Total length.
    Notes
    -----
    - This function is for a Network representing a single structural element, i.e. with two end-points (leaves).
    """

    start = rs.PointCoordinates(rs.ObjectsByLayer(layer)[0])

    return network_order(start=start, structure=structure, network=network)


def plot_reaction_forces(structure, step, layer=None, scale=1.0):
    """
    Plots reaction forces for the Structure analysis results.

    Parameters
    ----------
    structure : obj
        Structure object.
    step : str
        Name of the Step.
    layer : str
        Layer name for plotting.
    scale : float
        Scale of the arrows.

    Returns
    -------
    None

    """

    if not layer:
        layer = '{0}-{1}'.format(step, 'reactions')

    rs.CurrentLayer(rs.AddLayer(layer))
    rs.DeleteObjects(rs.ObjectsByLayer(layer))
    rs.EnableRedraw(False)

    rfx = structure.results[step]['nodal']['rfx']
    rfy = structure.results[step]['nodal']['rfy']
    rfz = structure.results[step]['nodal']['rfz']

    nkeys = rfx.keys()
    v = [scale_vector([rfx[i], rfy[i], rfz[i]], -scale * 0.001) for i in nkeys]
    rm = [length_vector(i) for i in v]
    rmax = max(rm)
    nodes = structure.nodes_xyz(nkeys)

    for i in nkeys:

        if rm[i] > 0.001:
            line = rs.AddLine(nodes[i], add_vectors(nodes[i], v[i]))
            rs.CurveArrows(line, 1)
            col = [int(j) for j in colorbar(rm[i] / rmax, input='float', type=255)]
            rs.ObjectColor(line, col)
            vector = [rfx[i], rfy[i], rfz[i]]
            name = json.dumps({'rfx': rfx[i], 'rfy': rfy[i], 'rfz': rfz[i], 'rfm': length_vector(vector)})
            rs.ObjectName(line, '_' + name)

    rs.CurrentLayer(rs.AddLayer('Default'))
    rs.LayerVisible(layer, False)
    rs.EnableRedraw(True)


def plot_concentrated_forces(structure, step, layer=None, scale=1.0):
    """
    Plots concentrated forces forces for the Structure analysis results.

    Parameters
    ----------
    structure : obj
        Structure object.
    step : str
        Name of the Step.
    layer : str
        Layer name for plotting.
    scale : float
        Scale of the arrows.

    Returns
    -------
    None

    """

    if not layer:
        layer = '{0}-{1}'.format(step, 'forces')
    rs.CurrentLayer(rs.AddLayer(layer))
    rs.DeleteObjects(rs.ObjectsByLayer(layer))
    rs.EnableRedraw(False)

    cfx = structure.results[step]['nodal']['cfx']
    cfy = structure.results[step]['nodal']['cfy']
    cfz = structure.results[step]['nodal']['cfz']

    nkeys = cfx.keys()
    v = [scale_vector([cfx[i], cfy[i], cfz[i]], -scale * 0.001) for i in nkeys]
    rm = [length_vector(i) for i in v]
    rmax = max(rm)
    nodes = structure.nodes_xyz(nkeys)

    for i in nkeys:

        if rm[i]:
            line = rs.AddLine(nodes[i], add_vectors(nodes[i], v[i]))
            rs.CurveArrows(line, 1)
            col = [int(j) for j in colorbar(rm[i] / rmax, input='float', type=255)]
            rs.ObjectColor(line, col)
            vector = [cfx[i], cfy[i], cfz[i]]
            name = json.dumps({'cfx': cfx[i], 'cfy': cfy[i], 'cfz': cfz[i], 'cfm': length_vector(vector)})
            rs.ObjectName(line, '_' + name)

    rs.CurrentLayer(rs.AddLayer('Default'))
    rs.LayerVisible(layer, False)
    rs.EnableRedraw(True)


def plot_mode_shapes(structure, step, layer=None, scale=1.0, radius=1):
    """
    Plots modal shapes from structure.results.

    Parameters
    ----------
    structure : obj
        Structure object.
    step : str
        Name of the Step.
    layer : str
        Each mode will be placed in a layer with this string prefix.
    scale : float
        Scale displacements for the deformed plot.
    radius : float
        Radius of the pipe visualisation meshes.

    Returns
    -------
    None

    """

    if not layer:
        layer = step + '_mode_'

    try:
        it = structure.results[step]['frequencies']
    except Exception:
        it = structure.results[step]['info']['description']

    if isinstance(it, list):
        for c, fk in enumerate(it, 1):
            layerk = layer + str(c)
            plot_data(structure=structure, step=step, field='um', layer=layerk, scale=scale, mode=c, radius=radius, source=None)

    elif isinstance(it, dict):
        for mode, value in it.items():
            #print(mode, value)
            layerk = layer + str(mode)
            plot_data(structure=structure, step=step, field='um', layer=layerk, scale=scale, mode=mode, radius=radius, source=None)


def plot_volmesh(volmesh, layer=None, draw_cells=True):
    """
    Plot a volmesh datastructure.

    Parameters
    ----------
    volmesh : obj
        volmesh datastructure object.
    layer : str
        Layer name to draw on.
    draw_cells : bool
        Draw cells.

    Returns
    -------
    None

    """

    if layer:
        rs.CurrentLayer(layer)

    vkeys = sorted(list(volmesh.vertices()), key=int)
    vertices = [volmesh.vertex_coordinates(vkey) for vkey in vkeys]

    if draw_cells:
        meshes = []
        for ckey in volmesh.cell:
            faces = [volmesh.halfface_vertices(fk, ordered=True) for fk in volmesh.cell_halffaces(ckey)]
            meshes.append(rs.AddMesh(vertices, faces))
        return meshes

    else:
        faces = []
        for fk in volmesh.halfface:
            face = volmesh.halfface_vertices(fk, ordered=True)
            faces.append(face)
        mesh = rs.AddMesh(vertices, faces)
        return mesh


def plot_axes(xyz, e11, e22, e33, layer, sc=1):
    """
    Plots a set of axes.

    Parameters
    ----------
    xyz : list
        Origin of the axes.
    e11 : list
        Normalised first axis component [x1, y1, z1].
    e22 : list
        Normalised second axis component [x2, y2, z2].
    e33 : list
        Normalised third axis component [x3, y3, z3].
    layer : str
        Layer to plot on.
    sc : float
         Size of the axis lines.

    Returns
    -------
    None

    """

    ex = rs.AddLine(xyz, add_vectors(xyz, scale_vector(e11, sc)))
    ey = rs.AddLine(xyz, add_vectors(xyz, scale_vector(e22, sc)))
    ez = rs.AddLine(xyz, add_vectors(xyz, scale_vector(e33, sc)))

    rs.ObjectColor(ex, [255, 0, 0])
    rs.ObjectColor(ey, [0, 255, 0])
    rs.ObjectColor(ez, [0, 0, 255])
    rs.ObjectLayer(ex, layer)
    rs.ObjectLayer(ey, layer)
    rs.ObjectLayer(ez, layer)


def plot_data(structure, lstep, field='um', layer=None, scale=1.0, radius=0.05, cbar=[None, None], iptype='mean',
              nodal='mean', mode='', cbar_size=1, source=None):
    """
    Plots analysis results on the deformed shape of the Structure.

    Parameters
    ----------
    structure : obj
        Structure object.
    lstep, step : str
        Name of the time Step.
    field : str
        Field to plot, e.g. 'um', 'sxx', 'sm1'.
    layer : str
        Layer name for plotting.
    scale : float
        Scale on displacements for the deformed plot.
    radius : float
        Radius of the pipe visualisation meshes.
    cbar : list
        Minimum and maximum limits on the colorbar.
    iptype : str
        'mean', 'max' or 'min' of an element's integration point data.
    nodal : str
        'mean', 'max' or 'min' for nodal values.
    mode : int
        Mode or frequency number to plot, for modal, harmonic or buckling analysis.
    cbar_size : float
        Scale on the size of the colorbar.
    source : float
        description of the used model

    Returns
    -------
    None

    Notes
    -----
    - Pipe visualisation of line elements is not based on the element section.

    """
   
    step=lstep
    tic=time()
    if field in ['smaxp', 'smises']:
        nodal = 'max'
        iptype = 'max'

    elif field in ['sminp']:
        nodal = 'min'
        iptype = 'min'

    # Create and clear Rhino layer

    if not layer:
        if source == 'linel':
            if field=='ux':
                name='(global x-displacements, source=linel)'            
            elif field == 'uy':
                name='(global y-displacements, source=linel)'
            elif field == 'uz':
                name='(global z-displacements, source=linel)'
            elif field == 'um':
                name='(global magnitue-displacements, source=linel)'
            elif field == 'sf1':
                name='(membrane forces in local x-direction, source=linel)'
            elif field == 'sf2':
                name='(membrane forces in local y-direction, source=linel)'
            elif field == 'sf3':
                name='(shear forces in local xy-direction, source=linel)'
            elif field == 'sf4':
                name='(transverse shear forces on local x-plane, source=linel)'
            elif field == 'sf5':
                name='(transverse shear forces on local y-plane, source=linel)'
            elif field == 'sm1':
                name='(bending moments around local y-direction, source=linel)'
            elif field == 'sm2':
                name='(bending moments around local x-direction, source=linel)'
            elif field == 'sm3':
                name='(twisting moments in x- and y-directions, source=linel)'   
            
        if source == 'SMM':
            if field == 'as_xi_bot':
                name='(amount of reinf. at bottom cover (local z-value is positiv) in local xi-direction [mm^2/m], source=SMM)'              
            elif field == 'as_xi_top':
                name='(amount of reinf. at top cover (local z-value is negativ) in local xi-direction [mm^2/m], source=SMM)'    
            elif field == 'as_eta_bot':
                name='(amount of reinf. at bottom cover (local z-value is positiv) in local eta-direction [mm^2/m], source=SMM)'    
            elif field == 'as_eta_top':
                name='(amount of reinf. at top cover (local z-value is negativ) in local eta-direction [mm^2/m], source=SMM)'    
            elif field == 'as_z':
                name='(amount of reinf. at in local z-value [mm^2/m], source=SMM)'    
            elif field == 'CC_bot':
                name='(values for concrete failure bottom cover (local z-value is positiv): 0=no failure, 1=failure but in iteration, 2=failure, source=SMM)'    
            elif field == 'CC_top':
                name='(values for concrete failure top cover (local z-value is negativ): 0=no failure, 1=failure but in iteration, 2=failure, source=SMM)'    
            elif field == 'k_bot':
                name='(k-factor of the sandwichmodel bottom cover (local z-value is positiv), source=SMM)'    
            elif field == 'k_top':
                name='(k-factor of the sandwichmodel top cover (local z-value is negativ), source=SMM)'    
            elif field == 't_bot':
                name='(thickness of the sandwich cover bottom cover (local z-value is positiv), source=SMM)'    
            elif field == 't_top':
                name='(thickness of the sandwich cover top cover (local z-direction is negativ), source=SMM)' 
            elif field == 'psi_bot':
                name='(angle between reinforcement directions xi and eta bottom cover (local z-value is positiv), source=SMM)'    
            elif field == 'psi_top':
                name='(angle between reinforcement directions xi and eta top cover (local z-value is negativ), source=SMM)'    
            elif field == 'Fall_bot':
                name='(case sandwichmodel cover bottom cover (local z-value is positiv), source=SMM)'    
            elif field == 'Fall_top':
                name='(case sandwichmodel cover top cover (local z-value is negativ), source=SMM)'    
            elif field == 'm_cc_bot':
                name='(degree of utilization for concrete failure bottom cover (local z-value is positiv), source=SMM)'    
            elif field == 'm_cc_top':
                name='(degree of utilization for concrete failure top cover (local z-value is negativ), source=SMM)'    
            elif field == 'm_shear_c':
                name='(degree of utilization for concrete core, source=SMM)'        
            elif field == 'm_c_total':
                name='(max degree of utilization for concrete, source=SMM)'                              

        if source == 'CMMUsermat':
            if field=='ux':
                name='(global x-displacements, source=CMMUsermat)'            
            elif field == 'uy':
                name='(global y-displacements, source=CMMUsermat)'
            elif field == 'uz':
                name='(global z-displacements, source=CMMUsermat)'
            elif field == 'um':
                name='(global magnitue-displacements, source=CMMUsermat)'
            elif field == 'sf1':
                name='(membrane forces in local x-direction, source=CMMUsermat)'
            elif field == 'sf2':
                name='(membrane forces in local y-direction, source=CMMUsermat)'
            elif field == 'sf3':
                name='(shear forces in local xy-direction, source=CMMUsermat)'
            elif field == 'sf4':
                name='(transverse shear forces on local x-plane, source=CMMUsermat)'
            elif field == 'sf5':
                name='(transverse shear forces on local y-plane, source=CMMUsermat)'
            elif field == 'sm1':
                name='(bending moments around local y-direction, source=CMMUsermat)'
            elif field == 'sm2':
                name='(bending moments around local x-direction, source=CMMUsermat)'
            elif field == 'sm3':
                name='(twisting moments in x- and y-directions, source=CMMUsermat)' 

        layer = '{0}-{1}-{3}{2}'.format(step, field, mode,name)

    rs.CurrentLayer(rs.AddLayer(layer))
    rs.DeleteObjects(rs.ObjectsByLayer(layer))
    rs.EnableRedraw(False)

    # Node and element data

    nodes = structure.nodes_xyz()
    elements = [structure.elements[i].nodes for i in sorted(structure.elements, key=int)]
    
    nodal_data = structure.results[step]['nodal']
    nkeys = sorted(structure.nodes, key=int)

    ux = [nodal_data['ux{0}'.format(mode)][i] for i in nkeys]
    uy = [nodal_data['uy{0}'.format(mode)][i] for i in nkeys]
    uz = [nodal_data['uz{0}'.format(mode)][i] for i in nkeys]

    try:
        data = [nodal_data['{0}{1}'.format(field, mode)][i] for i in nkeys]
        dtype = 'nodal'        

    except(Exception):
        data = structure.results[step]['element'][field]
        #print(data)
        dtype = 'element'

    

    # Postprocess
    
    result = functions.postprocess(nodes, elements, ux, uy, uz, data, dtype, scale, cbar, 255, iptype, nodal)
    
    try:
        toc, U, cnodes, fabs, fscaled, celements, eabs = result
        #print('\n***** Data processed : {0} s *****'.format(toc))

        # Plot meshes

        mesh_faces = []
        line_faces = [[0, 4, 5, 1], [1, 5, 6, 2], [2, 6, 7, 3], [3, 7, 4, 0]]
        block_faces = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]
        tet_faces = [[0, 2, 1, 1], [1, 2, 3, 3], [1, 3, 0, 0], [0, 3, 2, 2]]

        for element, nodes in enumerate(elements):
                        
            n = len(nodes)
            
            if n == 2:
                #if source != 'SMM':
                u, v = nodes
                sp, ep = U[u], U[v]
                plane = rs.PlaneFromNormal(sp, subtract_vectors(ep, sp))
                xa = plane.XAxis
                ya = plane.YAxis
                r = radius
                xa_pr = scale_vector(xa, +r)
                xa_mr = scale_vector(xa, -r)
                ya_pr = scale_vector(ya, +r)
                ya_mr = scale_vector(ya, -r)
                pts = [add_vectors(sp, xa_pr), add_vectors(sp, ya_pr),
                    add_vectors(sp, xa_mr), add_vectors(sp, ya_mr),
                    add_vectors(ep, xa_pr), add_vectors(ep, ya_pr),
                    add_vectors(ep, xa_mr), add_vectors(ep, ya_mr)]
                guid = rs.AddMesh(pts, line_faces)

                if dtype == 'element':
                    col1 = col2 = celements[element]
                    
                elif dtype == 'nodal':
                    col1 = cnodes[u]
                    col2 = cnodes[v]

                rs.MeshVertexColors(guid, [col1] * 4 + [col2] * 4)

            elif n == 3:

                mesh_faces.append(nodes + [nodes[-1]])

            elif n == 4:

                if structure.elements[element].__name__ in ['ShellElement', 'MembraneElement']:
                    mesh_faces.append(nodes)
                else:
                    for face in tet_faces:
                        mesh_faces.append([nodes[i] for i in face])

            elif n == 8:

                for block in block_faces:
                    mesh_faces.append([nodes[i] for i in block])

        if mesh_faces:
            guid = rs.AddMesh(U, mesh_faces)
            rs.MeshVertexColors(guid, cnodes)

        # Plot colorbar

        xr, yr, _ = structure.node_bounds()
        yran = yr[1] - yr[0] if yr[1] - yr[0] else 1
        s = yran * 0.1 * cbar_size
        xmin = xr[1] + 3 * s
        ymin = yr[0]

        xl = [xmin, xmin + s]
        yl = [ymin + i * s for i in range(11)]
        verts = [[xi, yi, 0] for xi in xl for yi in yl]
        faces = [[i, i + 1, i + 12, i + 11] for i in range(10)]
        id = rs.AddMesh(verts, faces)

        y = [i[1] for i in verts]
        yn = yran * cbar_size
        colors = [colorbar(2 * (yi - ymin - 0.5 * yn) / yn, input='float', type=255) for yi in y]
        rs.MeshVertexColors(id, colors)

        h = 0.008 * s

        for i in range(5):

            x0 = xmin + 1.2 * s
            yu = ymin + (5.8 + i) * s
            yl = ymin + (3.8 - i) * s
            vu = float(+max(eabs, fabs) * (i + 1) / 5.)
            vl = float(-max(eabs, fabs) * (i + 1) / 5.)
            rs.AddText('{0:.5g}'.format(vu), [x0, yu, 0], height=h)
            rs.AddText('{0:.5g}'.format(vl), [x0, yl, 0], height=h)

        rs.AddText('0', [x0, ymin + 4.8 * s, 0], height=h)
        rs.AddText('Step:{0}   Field:{1}'.format(step, field), [xmin, ymin + 12 * s, 0], height=h)

        if mode != '':
            try:
                freq = str(round(structure.results[step]['frequencies'][mode - 1], 3))
                rs.AddText('Mode:{0}   Freq:{1}Hz'.format(mode, freq), [xmin, ymin - 1.5 * s, 0], height=h)
            except Exception:
                pass

        # Return to Default layer

        rs.CurrentLayer(rs.AddLayer('Default'))
        rs.LayerVisible(layer, False)
        rs.EnableRedraw(True)


    except Exception:
        print('\n***** Error encountered during data processing or plotting *****')

    toc2 = time() - tic
    print('Plot {1} results in Rhino successful in {0:.3f} s'.format(toc2,field))



def plot_principal_stresses(structure, step, shell_layer, cbar_size=1, scale=10, layer=None, numeric='no', values='all'):
    """
    Plots the principal stresses of the elements.

    Parameters
    ----------
    structure : obj
        Structure object.
    step : str
        Name of the Step.
    sp : str
        'sp1' or 'sp5' for stection point 1 or 5.
    stype : str
        'max' or 'min' for maximum or minimum principal stresses.
    scale : float
        Scale on the length of the line markers (usually 10^6).
    layer : str
        Layer name for plotting.

    Returns
    -------
    None

    Notes
    -----
    - Centroids are taken on the undeformed geometry.

    """
    # Einlesen der Daten
    data = structure.results[step]['GP'] 
    

    # Pro Element jeweils 4 mal (=Anzahl GP) abfullen   
    
    if shell_layer == 'top':
        fcc_eff=data['fcc_eff_top']     
        coor_intp_layer_x=data['coor_intp_layer_x_top']     
        coor_intp_layer_y=data['coor_intp_layer_y_top']     
        coor_intp_layer_z=data['coor_intp_layer_z_top']
        
        
    elif shell_layer == 'bot':  
        fcc_eff=data['fcc_eff_bot']     
        coor_intp_layer_x=data['coor_intp_layer_x_bot']     
        coor_intp_layer_y=data['coor_intp_layer_y_bot']     
        coor_intp_layer_z=data['coor_intp_layer_z_bot']               
        
    elem_nr=data['elem_nr_bot']
    loc_x_glob_x=data['loc_x_glob_x']        
    loc_x_glob_y=data['loc_x_glob_y'] 
    loc_x_glob_z=data['loc_x_glob_z'] 
    loc_y_glob_x=data['loc_y_glob_x']        
    loc_y_glob_y=data['loc_y_glob_y'] 
    loc_y_glob_z=data['loc_y_glob_z']  
    elem_typ=data['elem_typ']   
    k_riss=data['k_riss']  
    fcc=data['fcc']





    # Aufbau Layer
    # --------------------------------------------------------------------------
    x_3_loc=[]
    y_3_loc=[]
    x_1_loc=[]
    y_1_loc=[]



    # Berechnung der Hauptspannungen und dessen Richtungen (in lokalen Koordinaten)
    # --------------------------------------------------------------------------
    
    ew_top, ev_top, ew_bot, ev_bot, length_stress=functions.principal_stresses(data,kind='sigma')   # ew = Eigenwerte (Hauptspannungen), ev=eigenvektoren (Hauptspannungsrichtungen),
  
    if shell_layer == 'top':
        ew=ew_top
        ev=ev_top
    elif shell_layer == 'bot':  
        ew=ew_bot
        ev=ev_bot 

    min_stress=min(min(ew))
    max_stress=max(max(ew))

        
    # Ploten der Hauptspannung 3
    # --------------------------------------------------------------------------
    if values=='all' or values=='3':
            
        # Layer in Rhino erzeugen
        if not layer:
                layer = '{0}_principal_stresses_{1}_layer_direction_3'.format(step, shell_layer)
        rs.CurrentLayer(rs.AddLayer(layer))
        rs.DeleteObjects(rs.ObjectsByLayer(layer))
        rs.EnableRedraw(False)


        # Plot Colorbar
        # --------------------------------------------------------------------------
        cbar_size=1

        xr, yr, _ = structure.node_bounds()
        yran = yr[1] - yr[0] if yr[1] - yr[0] else 1
        s = yran * 0.1 * cbar_size
        xmin = xr[1] + 3 * s
        ymin = yr[0]

        xl = [xmin, xmin + s]
        yl = [ymin + i * s for i in range(11)]
        h = 0.006 * s

        for i in range(3):
            x0 = xmin + 1.2 * s
            yl = ymin + (3.8 - i) * s
            if i==0:
                id=rs.AddText('0 < s', [x0, yl, 0], height=h)            
                col=[255,0,0]
                rs.ObjectColor(id, col)
            elif i==1:
                id=rs.AddText('0 > s > f_ceff', [x0, yl, 0], height=h)
                col=[0,0,255]
                rs.ObjectColor(id, col)
            elif i==2:
                id=rs.AddText('fceff < s', [x0, yl, 0], height=h)            
                col=[0,0,0]
                rs.ObjectColor(id, col)        

        for b in range(length_stress): 
                
            elem_typ_GP=elem_typ[b]

            if elem_typ_GP == 1: # 1=shell 0=MPR or others
                coor_intp_layer_x_GP=coor_intp_layer_x[b]
                coor_intp_layer_y_GP=coor_intp_layer_y[b]
                coor_intp_layer_z_GP=coor_intp_layer_z[b]
                loc_x_glob_x_GP=loc_x_glob_x[b]
                loc_x_glob_y_GP=loc_x_glob_y[b]
                loc_x_glob_z_GP=loc_x_glob_z[b]
                loc_y_glob_x_GP=loc_y_glob_x[b]      
                loc_y_glob_y_GP=loc_y_glob_y[b] 
                loc_y_glob_z_GP=loc_y_glob_z[b]          

                f_cc_eff_GP=fcc_eff[b]    
                fcc_GP=fcc[b]
                # Bestimmung aktuellen Hauptspannung 3
                # TODO: axes anpassen; anpassung fur hauptspannung
    


                f2 = Frame([coor_intp_layer_x_GP, coor_intp_layer_y_GP, coor_intp_layer_z_GP], [loc_x_glob_x_GP, loc_x_glob_y_GP, loc_x_glob_z_GP], [loc_y_glob_x_GP, loc_y_glob_y_GP, loc_y_glob_z_GP])    
                T = Transformation.from_frame(f2)  
                ev_GP=ev[b]
                ew_GP=ew[b]       
                
                sig_c3_GP=min(ew_GP)
                index_c3_GP=ew_GP.index(sig_c3_GP)     
                
                # Vectors local coor
                x_3_loc=ev_GP[0][index_c3_GP]*sig_c3_GP*scale
                y_3_loc=ev_GP[1][index_c3_GP]*sig_c3_GP*scale
                v_plus = Vector(x_3_loc*0.5, y_3_loc*0.5,0.).transformed(T)
                v_minus = Vector(-x_3_loc*0.5, -y_3_loc*0.5,0.).transformed(T)
                
                # Calculation mechanical values
                fct=k_riss[b]*0.3*fcc[b]**(2/3)
                
                centroid=[coor_intp_layer_x_GP,coor_intp_layer_y_GP,coor_intp_layer_z_GP]
                
                # Plot lines (only if abs(sig_c3_GP) > 0.01)
                if sig_c3_GP >= 0.01 and sig_c3_GP <= fct: # Zug ungerissen 
                    id3 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col3=[255,0,0] # rot
                    rs.ObjectColor(id3, col3)
                    if numeric == 'yes':                               
                        id3_text=rs.AddTextDot(str(round(sig_c3_GP,1))+' '+'('+str(round(fct/sig_c3_GP,1))+')', centroid)
                        rs.ObjectColor(id3_text, col3) 
                    else:
                        pass     
                elif sig_c3_GP > fct: # Zug gerissen
                    #id3 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    #col3=[255,120,0] # Orange
                    #rs.ObjectColor(id3, col3)
                    #if numeric == 'yes':                               
                    #    id3_text=rs.AddTextDot(str(round(sig_c3_GP,1))+' '+'('+str(round(fct/sig_c3_GP,1))+')', centroid)
                    #    rs.ObjectColor(id3_text, col3) 
                    #else:
                    #    pass 
                    pass # due to a brittle stress-strain relationshipt under tension -> sig_c1_GP=0 and therefore v_minus and v_plus are zero (not plotable)
                elif sig_c3_GP < -0.01 and sig_c3_GP > -1*f_cc_eff_GP: # kein Bruch                                       
                    id3 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col3=[0,0,255]
                    rs.ObjectColor(id3, col3)
                    if numeric == 'yes':                               
                        id3_text=rs.AddTextDot(str(round(sig_c3_GP,1))+' '+'('+str(round(-1*f_cc_eff_GP/sig_c3_GP,1))+')', centroid)                        
                        rs.ObjectColor(id3_text, col3) 
                    else:
                        pass                          
                elif sig_c3_GP <= -1*f_cc_eff_GP: # Bruch                    
                    id3 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col3=[0,0,0]
                    rs.ObjectColor(id3, col3)
                    if numeric == 'yes':                        
                        id3_text=rs.AddTextDot(str(round(sig_c3_GP,1))+' '+'('+str(round(-1*f_cc_eff_GP/sig_c3_GP,1))+')', centroid)
                        rs.ObjectColor(id3_text, col3) 
                    else:
                        pass  
                else: 
                    pass
                
                    # Folnder Command in Rhion ausfuhren falls liniendicke nicht angezeigt wird _PrintDisplay _State=_Toggle _Enter
            else:
                pass
    
    
    # Ploten der Hauptspannung 1
    # --------------------------------------------------------------------------
    if values=='all' or values=='1':
            
        # Layer in Rhino erzeugen
        if not layer:
                layer = '{0}_principal_stresses_{1}_layer_direction_1'.format(step, shell_layer)
        rs.CurrentLayer(rs.AddLayer(layer))
        rs.DeleteObjects(rs.ObjectsByLayer(layer))
        rs.EnableRedraw(False)



        # Plot Colorbar
        # --------------------------------------------------------------------------
        cbar_size=1

        xr, yr, _ = structure.node_bounds()
        yran = yr[1] - yr[0] if yr[1] - yr[0] else 1
        s = yran * 0.1 * cbar_size
        xmin = xr[1] + 3 * s
        ymin = yr[0]

        xl = [xmin, xmin + s]
        yl = [ymin + i * s for i in range(11)]
        h = 0.006 * s

        for i in range(3):
            x0 = xmin + 1.2 * s
            yl = ymin + (3.8 - i) * s
            if i==0:
                id=rs.AddText('0 < s', [x0, yl, 0], height=h)            
                col=[255,0,0]
                rs.ObjectColor(id, col)
            elif i==1:
                id=rs.AddText('0 > s > f_ceff', [x0, yl, 0], height=h)
                col=[0,0,255]
                rs.ObjectColor(id, col)
            elif i==2:
                id=rs.AddText('fceff < s', [x0, yl, 0], height=h)            
                col=[0,0,0]
                rs.ObjectColor(id, col)        
                        
        for b in range(length_stress): 
            
            elem_typ_GP=elem_typ[b]

            if elem_typ_GP == 1: # 1=shell 0=MPR or others
                coor_intp_layer_x_GP=coor_intp_layer_x[b]
                coor_intp_layer_y_GP=coor_intp_layer_y[b]
                coor_intp_layer_z_GP=coor_intp_layer_z[b]
                loc_x_glob_x_GP=loc_x_glob_x[b]
                loc_x_glob_y_GP=loc_x_glob_y[b]
                loc_x_glob_z_GP=loc_x_glob_z[b]
                loc_y_glob_x_GP=loc_y_glob_x[b]      
                loc_y_glob_y_GP=loc_y_glob_y[b] 
                loc_y_glob_z_GP=loc_y_glob_z[b]          

                f_cc_eff_GP=fcc_eff[b]           

                # Bestimmung aktuellen Hauptspannung 1
                # TODO: axes anpassen; anpassung fur hauptspannung

                f2 = Frame([coor_intp_layer_x_GP, coor_intp_layer_y_GP, coor_intp_layer_z_GP], [loc_x_glob_x_GP, loc_x_glob_y_GP, loc_x_glob_z_GP], [loc_y_glob_x_GP, loc_y_glob_y_GP, loc_y_glob_z_GP])    
                T = Transformation.from_frame(f2)  
                ev_GP=ev[b]
                ew_GP=ew[b]       

                sig_c1_GP=max(ew_GP)
                index_c1_GP=ew_GP.index(sig_c1_GP)     
                
                # Vectors local coor
                x_1_loc=ev_GP[0][index_c1_GP]*sig_c1_GP*scale
                y_1_loc=ev_GP[1][index_c1_GP]*sig_c1_GP*scale
                v_plus = Vector(x_1_loc*0.5, y_1_loc*0.5,0.).transformed(T)
                v_minus = Vector(-x_1_loc*0.5, -y_1_loc*0.5,0.).transformed(T)
                
                # Calculation mechanical values
                fct=k_riss[b]*0.3*fcc[b]**(2/3)
                
                centroid=[coor_intp_layer_x_GP,coor_intp_layer_y_GP,coor_intp_layer_z_GP]
                
                # Plot lines (only if abs(sig_c1_GP) > 0.01)
                if sig_c1_GP >= 0.01 and sig_c1_GP <= fct: # Zug ungerissen                    
                    id1 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col1=[255,0,0] # rot
                    rs.ObjectColor(id1, col1)
                    if numeric == 'yes':                               
                        id1_text=rs.AddTextDot(str(round(sig_c1_GP,1))+' '+'('+str(round(fct/sig_c1_GP,1))+')', centroid)
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass     
                elif sig_c1_GP > fct: # Zug gerissen
                    #id1 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    #col1=[255,120,0] # Orange
                    #rs.ObjectColor(id1, col1)
                    #if numeric == 'yes':                               
                    #    id1_text=rs.AddTextDot(str(round(sig_c1_GP,1))+' '+'('+str(round(fct/sig_c1_GP,1))+')', centroid)
                    #    rs.ObjectColor(id1_text, col1) 
                    #else:
                        #pass 
                    pass # due to a brittle stress-strain relationshipt under tension -> sig_c1_GP=0 and therefore v_minus and v_plus are zero (not plotable)
                elif sig_c1_GP < -0.01 and sig_c1_GP > -1*f_cc_eff_GP: # kein Bruch                                                      
                    id1 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col1=[0,0,255]
                    rs.ObjectColor(id1, col1)
                    if numeric == 'yes':                               
                        id1_text=rs.AddTextDot(str(round(sig_c1_GP,1))+' '+'('+str(round(-1*f_cc_eff_GP/sig_c1_GP,1))+')', centroid)
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass                          
                elif sig_c1_GP <= -1*f_cc_eff_GP: # Bruch                    
                    id1 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col1=[0,0,0]
                    rs.ObjectColor(id1, col1)
                    if numeric == 'yes':                        
                        id1_text=rs.AddTextDot(str(round(sig_c1_GP,1))+' '+'('+str(round(-1*f_cc_eff_GP/sig_c1_GP,1))+')', centroid)
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass   
                else:
                    pass

            
                
                    # Folnder Command in Rhion ausfuhren falls liniendicke nicht angezeigt wird _PrintDisplay _State=_Toggle _Enter
            else:
                pass
    

def plot_principal_strains(structure, step, shell_layer, cbar_size=1, scale=1000, layer=None, numeric='no', values='all'):
    """
    Plots the principal strains of the elements.

    Parameters
    ----------
    structure : obj
        Structure object.
    step : str
        Name of the Step.
    sp : str
        'sp1' or 'sp5' for stection point 1 or 5.
    stype : str
        'max' or 'min' for maximum or minimum principal stresses.
    scale : float
        Scale on the length of the line markers (usually 10^6).
    layer : str
        Layer name for plotting.

    Returns
    -------
    None

    Notes
    -----
    - Centroids are taken on the undeformed geometry.

    """
    # Einlesen der Daten
    data = structure.results[step]['GP'] 
    # Pro Element jeweils 4 mal (=Anzahl GP) abfullen   
    
    if shell_layer == 'top':
        eps_bruch=data['eps_bruch']     
        coor_intp_layer_x=data['coor_intp_layer_x_top']     
        coor_intp_layer_y=data['coor_intp_layer_y_top']     
        coor_intp_layer_z=data['coor_intp_layer_z_top']
        
        
    elif shell_layer == 'bot':  
        eps_bruch=data['eps_bruch']     
        coor_intp_layer_x=data['coor_intp_layer_x_bot']     
        coor_intp_layer_y=data['coor_intp_layer_y_bot']     
        coor_intp_layer_z=data['coor_intp_layer_z_bot']               
   

    elem_nr=data['elem_nr_bot']
    loc_x_glob_x=data['loc_x_glob_x']        
    loc_x_glob_y=data['loc_x_glob_y'] 
    loc_x_glob_z=data['loc_x_glob_z'] 
    loc_y_glob_x=data['loc_y_glob_x']        
    loc_y_glob_y=data['loc_y_glob_y'] 
    loc_y_glob_z=data['loc_y_glob_z']  
    elem_typ=data['elem_typ']   
    ecu=data['ecu']   
    k_E=data['k_E']   
    k_riss=data['k_riss']  
    fcc=data['fcc']

    

    # Aufbau Layer
    # --------------------------------------------------------------------------
    x_3_loc=[]
    y_3_loc=[]
    x_1_loc=[]
    y_1_loc=[]

    # Berechnung der Hauptverzerrungen und dessen Richtungen (in lokalen Koordinaten)
    # --------------------------------------------------------------------------
    # test area

    ew_top, ev_top, ew_bot, ev_bot, length_=functions.principal_stresses(data,kind='eps')  # ew = Eigenwerte (Hauptspannungen), ev=eigenvektoren (Hauptspannungsrichtungen),
    
    
    if shell_layer == 'top':
        ew=ew_top
        ev=ev_top
    elif shell_layer == 'bot':  
        ew=ew_bot
        ev=ev_bot 

    min_strains=min(min(ew))
    max_strains=max(max(ew))

        
    # Ploten der Hauptverzerrung 3
    # --------------------------------------------------------------------------
    if values=='all' or values=='3':
            
        # Layer in Rhino erzeugen
        if not layer:
                layer = '{0}_principal_strains_{1}_layer_direction_3'.format(step, shell_layer)
        rs.CurrentLayer(rs.AddLayer(layer))
        rs.DeleteObjects(rs.ObjectsByLayer(layer))
        rs.EnableRedraw(False)


        # Plot Colorbar
        # --------------------------------------------------------------------------
        xr, yr, _ = structure.node_bounds()
        yran = yr[1] - yr[0] if yr[1] - yr[0] else 1
        s = yran * 0.1 * cbar_size
        xmin = xr[1] + 3 * s
        ymin = yr[0]

        xl = [xmin, xmin + s]
        yl = [ymin + i * s for i in range(11)]
        h = 0.006 * s

        for i in range(4):
            x0 = xmin + 1.2 * s
            yl = ymin + (3.8 - i) * s
            if i==0:
                id=rs.AddText('0 < eps < etu', [x0, yl, 0], height=h)            
                col=[255,0,0]
                rs.ObjectColor(id, col)
            elif i==1:
                id=rs.AddText('eps > etu', [x0, yl, 0], height=h)
                col=[255,120,0]
                rs.ObjectColor(id, col)
            elif i==2:            
                id=rs.AddText('0 > eps > ecu', [x0, yl, 0], height=h)            
                col=[0,0,255]
                rs.ObjectColor(id, col)
            elif i==3:            
                id=rs.AddText('ecu > eps', [x0, yl, 0], height=h)            
                col=[0,0,0]
                rs.ObjectColor(id, col)



        for b in range(length_): 
                       
                
            elem_typ_GP=elem_typ[b]
            

            if elem_typ_GP == 1: # 1=shell 0=MPR or others
                coor_intp_layer_x_GP=coor_intp_layer_x[b]
                coor_intp_layer_y_GP=coor_intp_layer_y[b]
                coor_intp_layer_z_GP=coor_intp_layer_z[b]
                loc_x_glob_x_GP=loc_x_glob_x[b]
                loc_x_glob_y_GP=loc_x_glob_y[b]
                loc_x_glob_z_GP=loc_x_glob_z[b]
                loc_y_glob_x_GP=loc_y_glob_x[b]      
                loc_y_glob_y_GP=loc_y_glob_y[b] 
                loc_y_glob_z_GP=loc_y_glob_z[b]          

                eps_bruch_GP=eps_bruch[b]           

                # Bestimmung aktuellen Hauptverzerrung 3

                f2 = Frame([coor_intp_layer_x_GP, coor_intp_layer_y_GP, coor_intp_layer_z_GP], [loc_x_glob_x_GP, loc_x_glob_y_GP, loc_x_glob_z_GP], [loc_y_glob_x_GP, loc_y_glob_y_GP, loc_y_glob_z_GP])    
                T = Transformation.from_frame(f2)  
                ev_GP=ev[b]
                ew_GP=ew[b]       
                
                eps_3_GP=min(ew_GP)
                index_3_GP=ew_GP.index(eps_3_GP)     
                
                # Vectors local coor
                x_3_loc=ev_GP[0][index_3_GP]*eps_3_GP*scale
                y_3_loc=ev_GP[1][index_3_GP]*eps_3_GP*scale
                v_plus = Vector(x_3_loc*0.5, y_3_loc*0.5,0.).transformed(T)
                v_minus = Vector(-x_3_loc*0.5, -y_3_loc*0.5,0.).transformed(T)
                
                # Calculation mechanical values
                E_c=k_E[b]*fcc[b]**(1/3)                
                fct=k_riss[b]*0.3*fcc[b]**(2/3)
                eps_ctu=fct/E_c
                centroid=[coor_intp_layer_x_GP,coor_intp_layer_y_GP,coor_intp_layer_z_GP]
                

                # Plot lines (only if abs(eps_3_GP) > 0.000000001)
                if eps_3_GP >= 0.000000001 and eps_3_GP <= eps_ctu: # Zug ungerissen
                    id3 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col3=[255,0,0] # rot
                    rs.ObjectColor(id3, col3)
                    if numeric == 'yes':                               
                        id3_text=rs.AddTextDot(str(round(eps_3_GP*1000,1))+' '+'('+str(round(eps_ctu/eps_3_GP,1))+')', centroid)
                        rs.ObjectColor(id3_text, col3) 
                    else:
                        pass     
                elif eps_3_GP > eps_ctu: # Zug gerissen
                    id3 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col3=[255,120,0] # Orange
                    rs.ObjectColor(id3, col3)
                    if numeric == 'yes':                               
                        id3_text=rs.AddTextDot(str(round(eps_3_GP*1000,1))+' '+'('+str(round(eps_ctu/eps_3_GP,1))+')', centroid)
                        rs.ObjectColor(id3_text, col3) 
                    else:
                        pass 
                elif eps_3_GP < -0.000000001 and eps_3_GP > ecu[b]: # kein Bruch                                       
                    id3 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col3=[0,0,255]
                    rs.ObjectColor(id3, col3)
                    if numeric == 'yes':                               
                        id3_text=rs.AddTextDot(str(round(eps_3_GP*1000,1))+' '+'('+str(round(ecu[b]/eps_3_GP,1))+')', centroid)
                        rs.ObjectColor(id3_text, col3) 
                    else:
                        pass                          
                elif eps_3_GP <= ecu[b]: # Bruch                    
                    id3 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col3=[0,0,0]
                    rs.ObjectColor(id3, col3)
                    if numeric == 'yes':                        
                        id3_text=rs.AddTextDot(str(round(eps_3_GP*1000,1))+' '+'('+str(round(ecu[b]/eps_3_GP,1))+')', centroid)
                        rs.ObjectColor(id3_text, col3) 
                    else:
                        pass   
                else:
                    pass                      
            else:
                pass
                

    # Ploten der Hauptverzerrung 1
    # --------------------------------------------------------------------------
    if values=='all' or values=='1':

        # Layer in Rhino erzeugen
        if not layer:
                layer = '{0}_principal_strains_{1}_layer_direction_1'.format(step, shell_layer)
        rs.CurrentLayer(rs.AddLayer(layer))
        rs.DeleteObjects(rs.ObjectsByLayer(layer))
        rs.EnableRedraw(False)

    
        # Plot Colorbar
        # --------------------------------------------------------------------------
        xr, yr, _ = structure.node_bounds()
        yran = yr[1] - yr[0] if yr[1] - yr[0] else 1
        s = yran * 0.1 * cbar_size
        xmin = xr[1] + 3 * s
        ymin = yr[0]

        xl = [xmin, xmin + s]
        yl = [ymin + i * s for i in range(11)]
        h = 0.006 * s

        for i in range(4):
            x0 = xmin + 1.2 * s
            yl = ymin + (3.8 - i) * s
            if i==0:
                id=rs.AddText('0 < eps < etu', [x0, yl, 0], height=h)            
                col=[255,0,0]
                rs.ObjectColor(id, col)
            elif i==1:
                id=rs.AddText('eps > etu', [x0, yl, 0], height=h)
                col=[255,120,0]
                rs.ObjectColor(id, col)
            elif i==2:            
                id=rs.AddText('0 > eps > ecu', [x0, yl, 0], height=h)            
                col=[0,0,255]
                rs.ObjectColor(id, col)
            elif i==3:            
                id=rs.AddText('ecu > eps', [x0, yl, 0], height=h)            
                col=[0,0,0]
                rs.ObjectColor(id, col)        

        for b in range(length_): 
            
            elem_typ_GP=elem_typ[b]

            if elem_typ_GP == 1: # 1=shell 0=MPR or others
                coor_intp_layer_x_GP=coor_intp_layer_x[b]
                coor_intp_layer_y_GP=coor_intp_layer_y[b]
                coor_intp_layer_z_GP=coor_intp_layer_z[b]
                loc_x_glob_x_GP=loc_x_glob_x[b]
                loc_x_glob_y_GP=loc_x_glob_y[b]
                loc_x_glob_z_GP=loc_x_glob_z[b]
                loc_y_glob_x_GP=loc_y_glob_x[b]      
                loc_y_glob_y_GP=loc_y_glob_y[b] 
                loc_y_glob_z_GP=loc_y_glob_z[b]          

                eps_bruch_GP=eps_bruch[b]             

                # Bestimmung aktuellen Hauptverzerrung 1

                f2 = Frame([coor_intp_layer_x_GP, coor_intp_layer_y_GP, coor_intp_layer_z_GP], [loc_x_glob_x_GP, loc_x_glob_y_GP, loc_x_glob_z_GP], [loc_y_glob_x_GP, loc_y_glob_y_GP, loc_y_glob_z_GP])    
                T = Transformation.from_frame(f2)  
                ev_GP=ev[b]
                ew_GP=ew[b]       

                eps_1_GP=max(ew_GP)
                index_1_GP=ew_GP.index(eps_1_GP)     
                
                # Vectors local coor
                x_1_loc=ev_GP[0][index_1_GP]*eps_1_GP*scale
                y_1_loc=ev_GP[1][index_1_GP]*eps_1_GP*scale
                v_plus = Vector(x_1_loc*0.5, y_1_loc*0.5,0.).transformed(T)
                v_minus = Vector(-x_1_loc*0.5, -y_1_loc*0.5,0.).transformed(T)
                
                # Calculation mechanical values
                E_c=k_E[b]*fcc[b]**(1/3)
                fct=k_riss[b]*0.3*fcc[b]**(2/3)
                eps_ctu=fct/E_c
                centroid=[coor_intp_layer_x_GP,coor_intp_layer_y_GP,coor_intp_layer_z_GP]
                
                # Plot lines (only if abs(eps_3_GP) > 0.000000001)
                if eps_1_GP >= 0.000000001 and eps_1_GP <= eps_ctu: # Zug ungerissen
                    id1 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col1=[255,0,0] # rot
                    rs.ObjectColor(id1, col1)
                    if numeric == 'yes':                               
                        id1_text=rs.AddTextDot(str(round(eps_1_GP*1000,1))+' '+'('+str(round(eps_ctu/eps_1_GP,1))+')', centroid)
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass     
                elif eps_1_GP > eps_ctu: # Zug gerissen
                    id1 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col1=[255,120,0] # Orange
                    rs.ObjectColor(id1, col1)
                    if numeric == 'yes':                               
                        id1_text=rs.AddTextDot(str(round(eps_1_GP*1000,1))+' '+'('+str(round(eps_ctu/eps_1_GP,1))+')', centroid)
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass 
                elif eps_1_GP < -0.000000001 and eps_1_GP > ecu[b]: # kein Bruch                    
                    id1 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col1=[0,0,255]
                    rs.ObjectColor(id1, col1)
                    if numeric == 'yes':                               
                        id1_text=rs.AddTextDot(str(round(eps_1_GP*1000,1))+' '+'('+str(round(ecu[b]/eps_1_GP,1))+')', centroid)
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass                          
                elif eps_1_GP <= ecu[b]: # Bruch
                    id1 = rs.AddLine(add_vectors(centroid, v_minus), add_vectors(centroid, v_plus))   
                    col1=[0,0,0]
                    rs.ObjectColor(id1, col1)
                    if numeric == 'yes':                        
                        id1_text=rs.AddTextDot(str(round(eps_1_GP*1000,1))+' '+'('+str(round(ecu[b]/eps_1_GP,1))+')', centroid)
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass                         
                else:
                    pass
            else:
                pass



def plot_steel_stresses(structure, step, Reinf_layer, cbar_size=1, scale=1, layer=None, numeric='no'):
    """
    Plots the steel stresses on GP

    Parameters
    ----------
    structure : obj
        Structure object.
    step : str
        Name of the Step.
    sp : str
        'sp1' or 'sp5' for stection point 1 or 5.
    stype : str
        'max' or 'min' for maximum or minimum principal stresses.
    scale : float
        Scale on the length of the line markers (usually 10^6).
    layer : str
        Layer name for plotting.

    Returns
    -------
    None

    Notes
    -----
    - Centroids are taken on the undeformed geometry.

    """

    # Einlesen der Daten
    # --------------------------------------------------------------------------    
    data = structure.results[step]['GP'] 

        
        
    # Pro Element jeweils 4 mal (=Anzahl GP) abfullen   
    if Reinf_layer == 'RL_1':
        sig_sr_layer=data['sig_sr_1L']    
        coor_x_steel_layer=data['coor_x_sig_sr_1L']     
        coor_y_steel_layer=data['coor_y_sig_sr_1L']     
        coor_z_steel_layer=data['coor_z_sig_sr_1L']   
        psi=data['psi_1L']      
        fsy=data['fsy_1L']  
        fsu=data['fsu_1L']  

    elif Reinf_layer == 'RL_2':
        sig_sr_layer=data['sig_sr_2L']    
        coor_x_steel_layer=data['coor_x_sig_sr_2L']     
        coor_y_steel_layer=data['coor_y_sig_sr_2L']     
        coor_z_steel_layer=data['coor_z_sig_sr_2L']   
        psi=data['psi_2L']    
        fsy=data['fsy_2L']  
        fsu=data['fsu_2L']  

    elif Reinf_layer == 'RL_3':
        sig_sr_layer=data['sig_sr_3L']    
        coor_x_steel_layer=data['coor_x_sig_sr_3L']     
        coor_y_steel_layer=data['coor_y_sig_sr_3L']     
        coor_z_steel_layer=data['coor_z_sig_sr_3L']   
        psi=data['psi_3L'] 
        fsy=data['fsy_3L']  
        fsu=data['fsu_3L']      

    elif Reinf_layer == 'RL_4':
        sig_sr_layer=data['sig_sr_4L']    
        coor_x_steel_layer=data['coor_x_sig_sr_4L']     
        coor_y_steel_layer=data['coor_y_sig_sr_4L']     
        coor_z_steel_layer=data['coor_z_sig_sr_4L']   
        psi=data['psi_4L']    
        fsy=data['fsy_4L']  
        fsu=data['fsu_4L']  
        
    loc_x_glob_x=data['loc_x_glob_x']        
    loc_x_glob_y=data['loc_x_glob_y'] 
    loc_x_glob_z=data['loc_x_glob_z'] 
    loc_y_glob_x=data['loc_y_glob_x']        
    loc_y_glob_y=data['loc_y_glob_y'] 
    loc_y_glob_z=data['loc_y_glob_z']  
    elem_typ=data['elem_typ']   
    max_sig_sr=max(sig_sr_layer)    

    # Aufbau Layer
    # --------------------------------------------------------------------------
    if not layer:
            layer = '{0}_steel_stresses_{1}_Reinf_layer'.format(step, Reinf_layer)
    rs.CurrentLayer(rs.AddLayer(layer))
    rs.DeleteObjects(rs.ObjectsByLayer(layer))
    rs.EnableRedraw(False)


    # Plot Colorbar
    # --------------------------------------------------------------------------
    xr, yr, _ = structure.node_bounds()
    yran = yr[1] - yr[0] if yr[1] - yr[0] else 1
    s = yran * 0.1 * cbar_size
    xmin = xr[1] + 3 * s
    ymin = yr[0]

    xl = [xmin, xmin + s]
    yl = [ymin + i * s for i in range(11)]
    h = 0.006 * s

    for i in range(3):
        x0 = xmin + 1.2 * s
        yl = ymin + (3.8 - i) * s
        if i==0:
            id=rs.AddText('0 < s_sr < fsy', [x0, yl, 0], height=h)            
            col=[0,0,0]
            rs.ObjectColor(id, col)
        elif i==1:
            id=rs.AddText('fsy < s_sr < fsu', [x0, yl, 0], height=h)
            col=[255,165,0]
            rs.ObjectColor(id, col)
        elif i==2:            
            id=rs.AddText('fsu < s_sr', [x0, yl, 0], height=h)            
            col=[255,0,0]
            rs.ObjectColor(id, col) 


    # Ploten der Stahlspannungen
    # --------------------------------------------------------------------------
    length_steel_stress=len(sig_sr_layer) 

    for b in range(length_steel_stress): 
        
        elem_typ_GP=elem_typ[b]
        
        if elem_typ_GP == 1: # 1=shell 0=MPR or others
            sig_sr_GP=sig_sr_layer[b]
            
            if sig_sr_GP > 1: # plot only values higher 1 MPa 
                
                coor_intp_layer_x_GP=coor_x_steel_layer[b]
                coor_intp_layer_y_GP=coor_y_steel_layer[b]
                coor_intp_layer_z_GP=coor_z_steel_layer[b]
                angle=math.radians(psi[b])*-1 # Vorzeichen zur Berechnung der Darstellung gegenuber dem Usermat umkehren 
                loc_x_glob_x_GP=loc_x_glob_x[b]
                loc_x_glob_y_GP=loc_x_glob_y[b]
                loc_x_glob_z_GP=loc_x_glob_z[b]
                loc_y_glob_x_GP=loc_y_glob_x[b]      
                loc_y_glob_y_GP=loc_y_glob_y[b] 
                loc_y_glob_z_GP=loc_y_glob_z[b]          
                loc_x_glob_GP=[loc_x_glob_x_GP, loc_x_glob_y_GP, loc_x_glob_z_GP]
                loc_y_glob_GP=[loc_y_glob_x_GP, loc_y_glob_y_GP, loc_y_glob_z_GP]

                centroid=[coor_intp_layer_x_GP, coor_intp_layer_y_GP, coor_intp_layer_z_GP]  

                # step 1: Bestimmung der Vektoren aus den Stahlspannungen in globaler x-Richtung (Annahme: fikti alle Stahlspannung zeigen in x-Richtung)
                x_fiktiv=sig_sr_GP*scale

                y_fikiv=0      
                v_plus_fiktiv = Vector(x_fiktiv*0.5, y_fikiv*0.5,0.) 
                v_minus_fiktiv = Vector(-x_fiktiv*0.5, -y_fikiv*0.5,0.)                                                            
               
                # step 2: Transformation der globlane x-Richtung (world XY) in lokale Ebene (lokales Koordinantesystem)
                Frame_1 = Frame(centroid, loc_x_glob_GP, loc_y_glob_GP)  # Richtung lok koordinatensystem                
                T = Transformation.from_frame(Frame_1)                                  
                v_plus_fiktiv_trans=v_plus_fiktiv.transformed(T) # transformation from world XY to frame=lokales Koord.
                v_minus_fiktiv_trans=v_minus_fiktiv.transformed(T)  # transformation from world XY to frame=lokales Koord.
                
                # step 3: Normalvektor (lokale z-Achse in globlaen Koordianten) aus Kruezprodukt der lokalen x,y Achsen in globalen Koordianten bestimmen
                vector_cross=cross_vectors(loc_x_glob_GP,loc_y_glob_GP)
                vector_cross_normal=normalize_vector(vector_cross)

                # step 4: Rotation der nun transformieren Vektoren um den Bewehrungswinkle (angle) um die lokale z-Achse
                T = matrix_from_axis_and_angle(vector_cross_normal, angle) 
                v_plus_fiktiv_trans_rot=v_plus_fiktiv_trans.transformed(T) 
                v_minus_fiktiv_trans_rot=v_minus_fiktiv_trans.transformed(T) 

                # step 5: Ploten in Rhino (only if sig_sr_GP>0.01)
                if sig_sr_GP <= fsy[b] and sig_sr_GP > 0.01 :
                    id1 = rs.AddLine(add_vectors(centroid, v_minus_fiktiv_trans_rot), add_vectors(centroid, v_plus_fiktiv_trans_rot))                
                    col1=[0,0,0]
                    rs.ObjectColor(id1, col1)  
                    #rs.ObjectPrintWidth(id1, sig_sr_GP*0.01) # Dicke der Linie definieren (use _PrintDisplay in Rhino)                  
                    if numeric == 'yes':                        
                        id1_text=rs.AddTextDot(str(round(sig_sr_GP,1))+' '+'('+str(round(fsu[b]/sig_sr_GP,1))+')', centroid)
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass

                    #rs.ObjectPrintWidth(id1, sig_sr_GP*scale_line_width) # Dicke der Linie definieren (use _PrintDisplay in Rhino)
                elif sig_sr_GP > fsy[b] and sig_sr_GP <= fsu[b]:
                    id1 = rs.AddLine(add_vectors(centroid, v_minus_fiktiv_trans_rot), add_vectors(centroid, v_plus_fiktiv_trans_rot))                
                    col1=[255,165,0]
                    rs.ObjectColor(id1, col1)
                    #rs.ObjectPrintWidth(id1, sig_sr_GP*0.01) # Dicke der Linie definieren (use _PrintDisplay in Rhino) 
                    if numeric == 'yes':
                        id1_text=rs.AddTextDot(str(round(sig_sr_GP,1))+' '+'('+str(round(fsu[b]/sig_sr_GP,1))+')', centroid)                        
                    else:
                        pass              
                    #rs.ObjectPrintWidth(id1, sig_sr_GP*scale_line_width) # Dicke der Linie definieren (use _PrintDisplay in Rhino)
                elif sig_sr_GP > fsu[b]:
                    id1 = rs.AddLine(add_vectors(centroid, v_minus_fiktiv_trans_rot), add_vectors(centroid, v_plus_fiktiv_trans_rot))                
                    col1=[255,0,0]
                    rs.ObjectColor(id1, col1)
                    #rs.ObjectPrintWidth(id1, sig_sr_GP*0.01) # Dicke der Linie definieren (use _PrintDisplay in Rhino)
                    if numeric == 'yes':
                        id1_text=rs.AddTextDot(str(round(sig_sr_GP,1))+' '+'('+str(round(fsu[b]/sig_sr_GP,1))+')', centroid)                        
                    else:
                        pass                   
                    #rs.ObjectPrintWidth(id1, sig_sr_GP*scale_line_width) # Dicke der Linie definieren (use _PrintDisplay in Rhino)
                else:
                    pass

        else:
            pass



def plot_principal_shear(structure, step, field='shear', cbar_size=1, scale=10, layer=None, numeric='no', shear_verification='no', D_max=32, tau_cd=None):


    # extract data
    kmax = structure.element_count() # Anzahl Elemente, Startwert bei 1 nicht bei 0!
    data = structure.results[step]['element_info']
    
    loc_x_glob_x=data['elem_loc_x_glob_x'].values()        
    loc_x_glob_y=data['elem_loc_x_glob_y'].values()    
    loc_x_glob_z=data['elem_loc_x_glob_z'].values()    
    loc_y_glob_x=data['elem_loc_y_glob_x'].values()    
    loc_y_glob_y=data['elem_loc_y_glob_y'].values()    
    loc_y_glob_z=data['elem_loc_y_glob_z'].values()    
    
               
    # Extract data for SHEAR VERIFICATION
    if shear_verification=='yes':
        data_eps_x_06d_bot = structure.results[step]['GP']['eps_x_06d_bot']   
        data_eps_y_06d_bot = structure.results[step]['GP']['eps_y_06d_bot']   
        data_eps_xy_06d_bot = structure.results[step]['GP']['eps_xy_06d_bot']   
        
        data_eps_x_06d_top = structure.results[step]['GP']['eps_x_06d_top']   
        data_eps_y_06d_top = structure.results[step]['GP']['eps_y_06d_top']   
        data_eps_xy_06d_top = structure.results[step]['GP']['eps_xy_06d_top']   

        psi_1=structure.results[step]['element_info']['psi_1'].values() 
        psi_2=structure.results[step]['element_info']['psi_2'].values() 
        psi_3=structure.results[step]['element_info']['psi_3'].values() 
        psi_4=structure.results[step]['element_info']['psi_4'].values() 

        dm_1=structure.results[step]['element_info']['dm1'].values() 
        dm_2=structure.results[step]['element_info']['dm2'].values() 
        dm_3=structure.results[step]['element_info']['dm3'].values() 
        dm_4=structure.results[step]['element_info']['dm4'].values()   

        oo=structure.results[step]['element_info']['oo'].values() 
        uu=structure.results[step]['element_info']['uu'].values() 
        h_shell=structure.results[step]['element_info']['h_shell'].values() 
        LNr_d06_bot=structure.results[step]['element_info']['LNr_d06_bot'].values() 
        LNr_d06_top=structure.results[step]['element_info']['LNr_d06_top'].values()      
        fcc=structure.results[step]['element_info']['fcc'].values()      


    # Aufbau Layer
    # --------------------------------------------------------------------------
    if not layer:
            layer = '{0}_principial_shear'.format(step)
    rs.CurrentLayer(rs.AddLayer(layer))
    rs.DeleteObjects(rs.ObjectsByLayer(layer))
    rs.EnableRedraw(False)

    # Start 
    for k in range(kmax): 
        ele_type = structure.results[step]['element']['ele_type'][k].values()
        if ele_type[0] == 1.0:
            #sf4-->vx
            sf4 = structure.results[step]['element']['sf4'][k].values()
            vx = statistics.mean(list(sf4)) # Bildet den durchschnitt aller Integrationspunkte des Elements
            #sf5-->vy
            sf5 = structure.results[step]['element']['sf5'][k].values()
            vy = statistics.mean(list(sf5)) # Bildet den durchschnitt aller Integrationspunkte des Elements
            centroid= structure.element_centroid(element=k)   
             #vx+vy-->v0
            v0 = (vx**2+vy**2)**0.5
            
            angle_v0=math.atan(vy/vx) # Vorzeichen zur Berechnung der Darstellung gegenuber dem Usermat umkehren 

            loc_x_glob_GP=[loc_x_glob_x[k], loc_x_glob_y[k], loc_x_glob_z[k]]
            loc_y_glob_GP=[loc_y_glob_x[k], loc_y_glob_y[k], loc_y_glob_z[k]]          

            # step 1: Bestimmung der Vektoren aus den Stahlspannungen in globaler x-Richtung (Annahme: fikti alle Stahlspannung zeigen in x-Richtung)
            x_fiktiv=v0*scale
            y_fikiv=0      
            v_plus_fiktiv = Vector(x_fiktiv*0.5, y_fikiv*0.5,0.) 
            v_minus_fiktiv = Vector(-x_fiktiv*0.5, -y_fikiv*0.5,0.)                                                            
               
            # step 2: Transformation der globlane x-Richtung (world XY) in lokale Ebene (lokales Koordinantesystem)
            Frame_1 = Frame(centroid, loc_x_glob_GP, loc_y_glob_GP)  # Richtung lok koordinatensystem                
            T = Transformation.from_frame(Frame_1)                                  
            v_plus_fiktiv_trans=v_plus_fiktiv.transformed(T) # transformation from world XY to frame=lokales Koord.
            v_minus_fiktiv_trans=v_minus_fiktiv.transformed(T)  # transformation from world XY to frame=lokales Koord.
                
            # step 3: Normalvektor (lokale z-Achse in globlaen Koordianten) aus Kruezprodukt der lokalen x,y Achsen in globalen Koordianten bestimmen
            vector_cross=cross_vectors(loc_x_glob_GP,loc_y_glob_GP)
            vector_cross_normal=normalize_vector(vector_cross)

            # step 4: Rotation der nun transformieren Vektoren um den Bewehrungswinkle (angle) um die lokale z-Achse
            T = matrix_from_axis_and_angle(vector_cross_normal, angle_v0) 
            v_plus_fiktiv_trans_rot=v_plus_fiktiv_trans.transformed(T) 
            v_minus_fiktiv_trans_rot=v_minus_fiktiv_trans.transformed(T) 
            
            # step 5: Shear strength and verification (optional)
            if shear_verification=='yes':
                anlge_shear_ver=angle_v0 # Independon from a "main reinforcement"

                # Calculate mean values fro eps_06d bot and top (original values on GP and not at element midpoint!)    
                k_eps=k*4   
                eps_x_06d_bot_mean=(data_eps_x_06d_bot[k_eps]+data_eps_x_06d_bot[k_eps+1]+data_eps_x_06d_bot[k_eps+2]+data_eps_x_06d_bot[k_eps+3])/2
                eps_y_06d_bot_mean=(data_eps_y_06d_bot[k_eps]+data_eps_y_06d_bot[k_eps+1]+data_eps_y_06d_bot[k_eps+2]+data_eps_y_06d_bot[k_eps+3])/2
                eps_xy_06d_bot_mean=(data_eps_xy_06d_bot[k_eps]+data_eps_xy_06d_bot[k_eps+1]+data_eps_xy_06d_bot[k_eps+2]+data_eps_xy_06d_bot[k_eps+3])/2
                eps_06d_bot_mean=eps_x_06d_bot_mean*math.cos(anlge_shear_ver)**2+eps_y_06d_bot_mean*math.sin(anlge_shear_ver)**2+eps_xy_06d_bot_mean*math.cos(anlge_shear_ver)*math.sin(anlge_shear_ver)

                eps_x_06d_top_mean=(data_eps_x_06d_top[k_eps]+data_eps_x_06d_top[k_eps+1]+data_eps_x_06d_top[k_eps+2]+data_eps_x_06d_top[k_eps+3])/2
                eps_y_06d_top_mean=(data_eps_y_06d_top[k_eps]+data_eps_y_06d_top[k_eps+1]+data_eps_y_06d_top[k_eps+2]+data_eps_y_06d_top[k_eps+3])/2
                eps_xy_06d_top_mean=(data_eps_xy_06d_top[k_eps]+data_eps_xy_06d_top[k_eps+1]+data_eps_xy_06d_top[k_eps+2]+data_eps_xy_06d_top[k_eps+3])/2
                eps_06d_top_mean=eps_x_06d_top_mean*math.cos(anlge_shear_ver)**2+eps_y_06d_top_mean*math.sin(anlge_shear_ver)**2+eps_xy_06d_top_mean*math.cos(anlge_shear_ver)*math.sin(anlge_shear_ver)
                
                # Decision which values (top or bot) is decisive and calculation of the corresponding d_v (midpoint element)
                if eps_06d_top_mean > 0 and eps_06d_bot_mean < 0:
                    eps_06d_mean=eps_06d_top_mean # eps from top 
                    d_s1_top = h_shell[k]-uu[k]-dm_1[k]/2 
                    d_s2_top = h_shell[k]-uu[k]-dm_1[k]-dm_2[k]/2
                    d_v=(d_s1_top+d_s2_top)/2 # from top
                    Nr_layer=LNr_d06_top[k]
                elif eps_06d_top_mean > 0 and eps_06d_bot_mean > 0:
                    if eps_06d_top_mean > eps_06d_bot_mean:
                        eps_06d_mean=eps_06d_top_mean    
                        d_s1_top = h_shell[k]-uu[k]-dm_1[k]/2 
                        d_s2_top = h_shell[k]-uu[k]-dm_1[k]-dm_2[k]/2
                        d_v=(d_s1_top+d_s2_top)/2 # from top    
                        Nr_layer=LNr_d06_top[k]        
                    elif eps_06d_top_mean < eps_06d_bot_mean:
                        eps_06d_mean=eps_06d_bot_mean 
                        d_s4_bot = h_shell[k]-oo[k]-dm_4[k]/2 
                        d_s3_bot = h_shell[k]-oo[k]-dm_4[k]-dm_3[k]/2
                        d_v=(d_s4_bot+d_s3_bot)/2 # from bot 
                        Nr_layer=LNr_d06_bot[k]                   
                elif eps_06d_top_mean <= 0 and eps_06d_bot_mean <= 0:
                    eps_06d_mean=None #not defined as positiv
                elif eps_06d_top_mean < 0 and eps_06d_bot_mean > 0:
                    eps_06d_mean=eps_06d_bot_mean
                    d_s4_bot = h_shell[k]-oo[k]-dm_4[k]/2 
                    d_s3_bot = h_shell[k]-oo[k]-dm_4[k]-dm_3[k]/2
                    d_v=(d_s4_bot+d_s3_bot)/2 # from bot  
                    Nr_layer=LNr_d06_bot[k]

                # Calculation of shear resistance accorindg to SIA 269              
                if eps_06d_mean != None:    
                    k_g=48/(16+D_max)
                    k_d=1/(1+2.5*eps_06d_mean*d_v*k_g)                
                    v_rd=k_d*tau_cd*d_v                    
                else:
                    v_rd=None


            # step 6: Ploten in Rhino (only if sig_sr_GP>0.01)
            if shear_verification=='no':
                id = rs.AddLine(add_vectors(centroid, v_minus_fiktiv_trans_rot), add_vectors(centroid, v_plus_fiktiv_trans_rot))                
                col=[0,0,0]
                rs.ObjectColor(id, col)  
                if numeric == 'yes':                        
                    id1_text=rs.AddTextDot(str(round(v0,1)), centroid)
                    rs.ObjectColor(id1_text, col) 

            elif shear_verification=='yes':
                if v_rd != None: 
                    if v0 <= v_rd:                    
                        id = rs.AddLine(add_vectors(centroid, v_minus_fiktiv_trans_rot), add_vectors(centroid, v_plus_fiktiv_trans_rot))                
                        col=[0,0,0]
                        rs.ObjectColor(id, col)  
                    elif v0 > v_rd:
                        id = rs.AddLine(add_vectors(centroid, v_minus_fiktiv_trans_rot), add_vectors(centroid, v_plus_fiktiv_trans_rot))                
                        col=[255,0,0]
                        rs.ObjectColor(id, col)                                                              
                    if numeric == 'yes':                        
                        id1_text=rs.AddTextDot(str(round(v0,1))+' '+'('+str(round(v_rd/v0,1))+')', centroid)
                        col1=[0,0,0]
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass
                else:
                    id = rs.AddLine(add_vectors(centroid, v_minus_fiktiv_trans_rot), add_vectors(centroid, v_plus_fiktiv_trans_rot))                
                    col=[0,0,0]
                    rs.ObjectColor(id, col)   
                    if numeric == 'yes':                        
                        id1_text=rs.AddTextDot(str(round(v0,1))+' '+'(not def)', centroid)
                        col1=[0,0,0]
                        rs.ObjectColor(id1_text, col1) 
                    else:
                        pass                                      

            # Plot Colorbar
            
            # --------------------------------------------------------------------------
            xr, yr, _ = structure.node_bounds()
            yran = yr[1] - yr[0] if yr[1] - yr[0] else 1
            s = yran * 0.1 * cbar_size
            xmin = xr[1] + 3 * s
            ymin = yr[0]

            xl = [xmin, xmin + s]
            yl = [ymin + i * s for i in range(11)]
            h = 0.006 * s

            for i in range(2):
                x0 = xmin + 1.2 * s
                yl = ymin + (3.8 - i) * s
                if i==0:
                    id=rs.AddText('v0 <= vrd', [x0, yl, 0], height=h)            
                    col=[0,0,0]
                    rs.ObjectColor(id, col)
                elif i==1:
                    id=rs.AddText('v0 > vrd', [x0, yl, 0], height=h)
                    col=[255,0,0]
                    rs.ObjectColor(id, col)


def plot_voxels(structure, step, field='smises', cbar=[None, None], iptype='mean', nodal='mean', vdx=None, mode=''):
    """
    Voxel 4D visualisation.

    Parameters
    ----------
    structure : obj
        Structure object.
    step : str
        Name of the Step.
    field : str
        Field to plot, e.g. 'smises'.
    cbar : list
        Minimum and maximum limits on the colorbar.
    iptype : str
        'mean', 'max' or 'min' of an element's integration point data.
    nodal : str
        'mean', 'max' or 'min' for nodal values.
    vdx : float
        Voxel spacing.
    mode : int
        mode or frequency number to plot, in case of modal, harmonic or buckling analysis.

    Returns
    -------
    None

    """

    # Node and element data

    xyz = structure.nodes_xyz()
    elements = [structure.elements[i].nodes for i in sorted(structure.elements, key=int)]
    nodal_data = structure.results[step]['nodal']
    nkeys = sorted(structure.nodes, key=int)

    ux = [nodal_data['ux{0}'.format(mode)][i] for i in nkeys]
    uy = [nodal_data['uy{0}'.format(mode)][i] for i in nkeys]
    uz = [nodal_data['uz{0}'.format(mode)][i] for i in nkeys]

    try:
        data = [nodal_data[field + str(mode)][key] for key in nkeys]
        dtype = 'nodal'

    except(Exception):
        data = structure.results[step]['element'][field]
        dtype = 'element'

    # Postprocess

    result = functions.postprocess(xyz, elements, ux, uy, uz, data, dtype, 1, cbar, 255, iptype, nodal)

    try:
        toc, U, cnodes, fabs, fscaled, celements, eabs = result
        print('\n***** Data processed : {0} s *****'.format(toc))

    except Exception:
        print('\n***** Error post-processing *****')

    try:
        functions.plotvoxels(values=fscaled, U=U, vdx=vdx)
        print('\n***** Voxels finished *****')

    except Exception:
        print('\n***** Error plotting voxels *****')


def weld_meshes_from_layer(layer_input, layer_output):
    """
    Grab meshes on an input layer and weld them onto an output layer.

    Parameters
    ----------
    layer_input : str
        Layer containing the Rhino meshes to weld.
    layer_output : str
        Layer to plot single welded mesh.

    Returns
    -------
    None

    """

    print('Welding meshes on layer:{0}'.format(layer_input))

    mdl = Structure(path=' ')

    add_nodes_elements_from_layers(mdl, mesh_type='ShellElement', layers=layer_input)

    faces = []

    for element in mdl.elements.values():
        enodes = element.nodes

        if len(enodes) == 3:
            enodes.append(enodes[-1])

        if len(enodes) == 4:
            faces.append(enodes)

    rs.DeleteObjects(rs.ObjectsByLayer(layer_output))
    rs.CurrentLayer(layer_output)
    rs.AddMesh(mdl.nodes_xyz(), faces)

# ==============================================================================
# Debugging
# ==============================================================================


if __name__ == "__main__":

    pass
