"""
********************************************************************************
fea
********************************************************************************

.. currentmodule:: compas_fea.fea

The compas_fea package with struceng lib supports Ansys as analysis backends.


Classes
=======

.. autosummary::
    :toctree: generated/

    Writer


Backends
========


ansys
-----

.. currentmodule:: compas_fea.fea.ansys

.. autosummary::
    :toctree: generated/

    input_generate
    make_command_file_static
    make_command_file_modal
    make_command_file_harmonic
    ansys_launch_process
    ansys_launch_process_extract
    delete_result_files
    extract_rst_data
    write_results_from_rst
    load_to_results


"""
from __future__ import absolute_import

from .writer import Writer

# Author(s): Compas/Compas FEA Team, Marius  Weber (ETHZ, HSLU T&A)

__all__ = [
    'Writer'
]
