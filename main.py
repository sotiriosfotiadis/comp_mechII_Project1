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


def main():
    # 1. Διάβασμα του .feb / .semfe αρχείου
    feb_file = "validation.semfe"   # ή όποιο όνομα έχεις
    nodes, elems, materials, k, bcs = read_input_file(feb_file)

    # 2. Συναρμολόγηση global K και f
    K = assemble_global(nodes, elems, k)
    f = np.zeros(nodes.shape[0])

    # 3. Εφαρμογή Neumann (heat flux) BCs
    #    bcs['heat_flux'] είναι λίστα (elem_id, edge_id, q)
    if bcs.get("heat_flux"):
        f = apply_heat_flux(f, nodes, elems, bcs["heat_flux"])

    # 4. Εφαρμογή Convection (Robin) BCs
    #    bcs['convection'] είναι λίστα (elem_id, edge_id, h, Tinf)
    if bcs.get("convection"):
        K, f = apply_convection(K, f, nodes, elems, bcs["convection"])

    # 5. Εφαρμογή Dirichlet (Temperature) BCs
    #    bcs['temperature'] είναι λίστα (node, value)
    if bcs.get("temperature"):
        bc_nodes  = [node for node, val in bcs["temperature"]]
        bc_values = [val  for node, val in bcs["temperature"]]
        K, f = apply_dirichlet(K, f, bc_nodes, bc_values)

    # 6. Λύση συστήματος
    u = solve_system(K, f)

    # 7. Post-processing
    #    (αν θες, ξεσχολιάζεις και το plot_mesh / plot_mesh_interactive)
    plot_mesh(nodes, elems, filename="mesh.png")
    plot_mesh_interactive(nodes, elems)

    plot_temperature_field(nodes, elems, u, filename="temperature_field.png")
    export_temperature_csv(nodes, u)


if __name__ == "__main__":
    main()