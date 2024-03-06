"""

********************************************************************************
cad
********************************************************************************




.. currentmodule:: compas_fea.cad

The compas_fea package for struceng lib supports Rhino in the frontend.


Rhino
=====

.. autosummary::
    :toctree: generated/
    :nosignatures:

    add_element_set
    add_node_set
    add_nodes_elements_from_layers
    add_sets_from_layers
    add_tets_from_mesh
    discretise_mesh
    mesh_extrude
    network_from_lines
    ordered_network
    plot_reaction_forces
    plot_concentrated_forces
    plot_mode_shapes
    plot_volmesh
    plot_axes
    plot_data
    plot_principal
    plot_voxels
    weld_meshes_from_layer

"""

# Author(s): Compas/Compas FEA Team, Marius  Weber (ETHZ, HSLU T&A)

from __future__ import absolute_import

import compas

if compas.RHINO:
    from .rhino import (  # noqa: F401
        add_element_set,
        add_node_set,
        add_nodes_elements_from_layers,
        add_sets_from_layers,
        add_tets_from_mesh,
        discretise_mesh,
        mesh_extrude,
        network_from_lines,
        ordered_network,
        plot_reaction_forces,
        plot_concentrated_forces,
        plot_mode_shapes,
        plot_volmesh,
        plot_axes,
        plot_data,
        plot_principal,
        plot_voxels,
        weld_meshes_from_layer,
    )
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
        'plot_principal',
        'plot_voxels',
        'weld_meshes_from_layer',
    ]

else:
    __all__ = []
