#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The SEMFE Heat Transfer Solver
Computational Mechanics

Main Script
"""
import numpy as np
from PreProcessor import read_input_file
from Solver import assemble_global, apply_convection, apply_dirichlet
from Solver import apply_heat_flux, solve_system
from PostProcessor import plot_mesh, plot_mesh_interactive, plot_temperature_field
from PostProcessor import export_temperature_csv


# Import model info
nodes, elems, materials, k, bcs = read_input_file('Ex1.semfe')

# Check Mesh Quality
plot_mesh_interactive(nodes, elems, show=True, filename='interactive_mesh2.html')

# Assemble global
K = assemble_global(nodes, elems, k=k)

# Apply BCs
bc_nodes = [node for node, val in bcs['temperature']]
bc_values = [val for node, val in bcs['temperature']]
Kmod, fmod = apply_dirichlet(K, np.zeros(nodes.shape[0]), bc_nodes, bc_values)
fmod       = apply_heat_flux(fmod, nodes, elems, bc_values)
Kmod, fmod = apply_convection(K, fmod, nodes, elems, bc_values)

#
u = solve_system(Kmod, fmod)

# Call it in main
plot_temperature_field(nodes, elems, u, filename='temperature_field.png')
export_temperature_csv(nodes, u)
