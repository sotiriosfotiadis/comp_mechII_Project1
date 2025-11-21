# -*- coding: utf-8 -*-
"""
Chimney mesh + SEMFE XML creator
- structured grid
- inner rectangular hole
- triangular elements (tri3)
"""

from collections import defaultdict

# -----------------------------
# ΓΕΩΜΕΤΡΙΑ ΚΑΜΙΝΑΔΑΣ
# -----------------------------
OUTER_W = 0.8     # πλάτος εξωτερικού ορθογωνίου
OUTER_H = 0.6     # ύψος εξωτερικού ορθογωνίου

# οπή (εσωτερικό ορθογώνιο)
HOLE_X1, HOLE_X2 = 0.2, 0.6
HOLE_Y1, HOLE_Y2 = 0.2, 0.4

# ΥΛΙΚΟ
K_COND = 1.5      # θερμική αγωγιμότητα

# ΤΙΜΕΣ ΣΥΝΟΡΙΑΚΩΝ ΣΥΝΘΗΚΩΝ
T_INNER = 100.0   # θερμοκρασία στο εσωτερικό τοίχωμα
T_TOP   = 30.0    # θερμοκρασία στο πάνω εξωτερικό τοίχωμα
H_CONV  = 50.0    # συντελεστής συναγωγής στο δεξί εξωτερικό
T_INF   = 25.0    # θερμοκρασία ρευστού


# ---------------------------------------------------------
# 1. ΔΗΜΙΟΥΡΓΙΑ ΠΛΕΓΜΑΤΟΣ
# ---------------------------------------------------------
def generate_chimney_mesh(nx=20, ny=15):
    """
    Φτιάχνει structured ορθογώνιο grid (0..OUTER_W, 0..OUTER_H),
    πετάει τα cells που πέφτουν μέσα στο εσωτερικό ορθογώνιο (hole),
    και κάθε cell το κόβει σε 2 τρίγωνα.

    Επιστρέφει:
    - nodes: λίστα [(x,y), ...]
    - elems: λίστα [(n1,n2,n3), ...] με 0-based indices κόμβων
    """

    dx = OUTER_W / nx
    dy = OUTER_H / ny

    n_nodes_x = nx + 1
    n_nodes_y = ny + 1

    # ----- κόμβοι -----
    nodes = []
    for j in range(n_nodes_y):
        for i in range(n_nodes_x):
            x = i * dx
            y = j * dy
            nodes.append((x, y))

    # helper: index κόμβου (i,j) → global id
    def nid(i, j):
        return j * n_nodes_x + i

    # ----- στοιχεία (triangles) -----
    elems = []
    for j in range(ny):
        for i in range(nx):
            # κέντρο του cell
            xc = (i + 0.5) * dx
            yc = (j + 0.5) * dy

            # αν το κέντρο είναι μέσα στο εσωτερικό ορθογώνιο → skip
            if (HOLE_X1 < xc < HOLE_X2) and (HOLE_Y1 < yc < HOLE_Y2):
                continue

            # κόμβοι cell (τετράπλευρο)
            n1 = nid(i,   j)
            n2 = nid(i+1, j)
            n3 = nid(i+1, j+1)
            n4 = nid(i,   j+1)

            # 2 τρίγωνα (n1,n2,n4) και (n2,n3,n4)
            elems.append((n1, n2, n4))
            elems.append((n2, n3, n4))

    return nodes, elems


# ---------------------------------------------------------
# 2. ΕΥΡΕΣΗ ΣΥΝΟΡΙΑΚΩΝ ΑΚΜΩΝ
# ---------------------------------------------------------
def find_boundary_edges(nodes, elems):
    """
    Βρίσκει όλες τις ακμές που ανήκουν σε μόνο 1 τρίγωνο → boundary edges.

    Επιστρέφει λίστα (a, b, elem_id, local_edge_id)
    """
    edge_map = defaultdict(list)

    for e_idx, (n1, n2, n3) in enumerate(elems):
        # τοπικές ακμές με local ID 1–3
        edges = [
            (n1, n2, 1),
            (n2, n3, 2),
            (n3, n1, 3),
        ]
        for a, b, local in edges:
            key = (min(a, b), max(a, b))  # χωρίς κατεύθυνση
            edge_map[key].append((e_idx, local))

    boundary_edges = []
    for (a, b), owners in edge_map.items():
        if len(owners) == 1:  # μόνο ένα στοιχείο → σύνορο
            e_idx, local = owners[0]
            boundary_edges.append((a, b, e_idx, local))

    return boundary_edges


# ---------------------------------------------------------
# 3. ΤΑΞΙΝΟΜΗΣΗ ΣΥΝΟΡΩΝ
# ---------------------------------------------------------
def classify_edges(nodes, bnd):
    """
    Χωρίζει τις συνοριακές ακμές σε:
    - bottom (y=0)
    - right_ext (x=OUTER_W)
    - top_ext (y=OUTER_H)
    - left_ext (x=0)  [μονωμένο, δεν το χρειαζόμαστε στο XML]
    - inner  (ακμές γύρω από την τρύπα)
    """
    tol = 1e-10
    bottom = []
    right_ext = []
    top_ext = []
    left_ext = []
    inner = []

    for (a, b, e_idx, loc) in bnd:
        x1, y1 = nodes[a]
        x2, y2 = nodes[b]

        if abs(y1) < tol and abs(y2) < tol:
            bottom.append((a, b, e_idx, loc))
        elif abs(x1 - OUTER_W) < tol and abs(x2 - OUTER_W) < tol:
            right_ext.append((a, b, e_idx, loc))
        elif abs(y1 - OUTER_H) < tol and abs(y2 - OUTER_H) < tol:
            top_ext.append((a, b, e_idx, loc))
        elif abs(x1) < tol and abs(x2) < tol:
            left_ext.append((a, b, e_idx, loc))
        else:
            inner.append((a, b, e_idx, loc))

    return bottom, right_ext, top_ext, left_ext, inner


# ---------------------------------------------------------
# 4. ΚΑΘΑΡΙΣΜΟΣ / ΕΠΑΝΑΡΙΘΜΗΣΗ ΚΟΜΒΩΝ
# ---------------------------------------------------------
def cleanup(nodes, elems, bottom, right, top, left, inner):
    """
    1) Πετάει τυχόν κόμβους που δεν χρησιμοποιούνται σε κανένα στοιχείο.
    2) Κάνει reindex ώστε οι κόμβοι να είναι 0..N-1 χωρίς κενά.
    3) Ενημερώνει αντίστοιχα elems & boundary edges.
    """
    # ποιοι κόμβοι χρησιμοποιούνται
    used = set()
    for n1, n2, n3 in elems:
        used.update([n1, n2, n3])

    used = sorted(list(used))
    mapping = {old: new for new, old in enumerate(used)}

    # νέοι κόμβοι (μόνο οι used, με νέα αρίθμηση)
    new_nodes = [nodes[i] for i in used]

    # νέα στοιχεία με reindexed κόμβους
    new_elems = [(mapping[a], mapping[b], mapping[c]) for (a, b, c) in elems]

    # helper για reindex στα edges
    def remap(edges):
        return [(mapping[a], mapping[b], e_idx, local)
                for (a, b, e_idx, local) in edges]

    return (
        new_nodes,
        new_elems,
        remap(bottom),
        remap(right),
        remap(top),
        remap(left),
        remap(inner),
    )


# ---------------------------------------------------------
# 5. ΓΡΑΨΙΜΟ XML (chimney.semfe)
# ---------------------------------------------------------
def write_xml(filename, nodes, elems, bottom, right, top, inner):
    """
    Γράφει αρχείο SEMFE σε μορφή συμβατή με το PreProcessor.
    - Tri3 elements
    - Dirichlet σε inner & top
    - Heat flux q=0 στο bottom
    - Convection στο right
    """

    with open(filename, "w", encoding="ISO-8859-1") as f:

        # HEADER
        f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        f.write('<SEMFE_spec>\n')
        f.write('  <Module type="heat conduction"/>\n\n')

        # MATERIALS
        f.write('  <Materials>\n')
        f.write('    <Material id="1" name="Brick">\n')
        f.write(f'      <conductivity>{K_COND}</conductivity>\n')
        f.write('    </Material>\n')
        f.write('  </Materials>\n\n')

        # GEOMETRY: NODES
        f.write('  <Geometry>\n')
        f.write('    <Nodes>\n')
        for i, (x, y) in enumerate(nodes, start=1):
            f.write(f'      <node id="{i}" x="{x}" y="{y}" z="0.0"/>\n')
        f.write('    </Nodes>\n\n')

        # GEOMETRY: ELEMENTS
        f.write('    <Elements type="tri3" name="mesh">\n')
        for e_id, (a, b, c) in enumerate(elems, start=1):
            # +1 γιατί το SEMFE χρησιμοποιεί 1-based ids
            f.write(f'      <elem id="{e_id}">{a+1} {b+1} {c+1}</elem>\n')
        f.write('    </Elements>\n')
        f.write('  </Geometry>\n\n')

        # BOUNDARY CONDITIONS
        f.write('  <BoundaryConditions>\n')

        # ---- DIRICHLET ----
        f.write('    <Boundary>\n')

        # Εσωτερικό τοίχωμα T = T_INNER
        inner_nodes = sorted({a for a, b, _, _ in inner} |
                             {b for a, b, _, _ in inner})
        for n in inner_nodes:
            f.write(f'      <temperature node="{n+1}" value="{T_INNER}"/>\n')

        # Πάνω εξωτερικό T = T_TOP
        top_nodes = sorted({a for a, b, _, _ in top} |
                           {b for a, b, _, _ in top})
        for n in top_nodes:
            f.write(f'      <temperature node="{n+1}" value="{T_TOP}"/>\n')

        f.write('    </Boundary>\n\n')

        # ---- HEAT FLUX bottom = 0 ----
        f.write('    <HeatFlux>\n')
        for (a, b, e_idx, local) in bottom:
            f.write(f'      <flux elem="{e_idx+1}" edge="{local}" value="0.0"/>\n')
        f.write('    </HeatFlux>\n\n')

        # ---- CONVECTION στο δεξί τοίχωμα ----
        f.write('    <Convection>\n')
        for (a, b, e_idx, local) in right:
            f.write(
                f'      <conv elem="{e_idx+1}" edge="{local}" '
                f'h="{H_CONV}" Tinf="{T_INF}"/>\n'
            )
        f.write('    </Convection>\n')

        f.write('  </BoundaryConditions>\n\n')

        # STEP
        f.write('  <Step name="step1" type="steady-state">\n')
        f.write('    <HeatSource>\n')
        f.write('    </HeatSource>\n')
        f.write('  </Step>\n')
        f.write('</SEMFE_spec>\n')

    print("XML saved:", filename)


# ---------------------------------------------------------
# MAIN
# ---------------------------------------------------------
if __name__ == "__main__":
    # 1) Δημιουργία πλέγματος
    nodes, elems = generate_chimney_mesh(nx=20, ny=15)

    # 2) Εντοπισμός συνοριακών ακμών
    bnd = find_boundary_edges(nodes, elems)

    # 3) Ταξινόμηση (κάτω, δεξί, πάνω, αριστερό, εσωτερικό)
    bottom, right, top, left, inner = classify_edges(nodes, bnd)

    # 4) Καθάρισμα & reindex όλων
    nodes, elems, bottom, right, top, left, inner = cleanup(
        nodes, elems, bottom, right, top, left, inner
    )

    # 5) Γράψιμο XML
    write_xml("chimney.semfe", nodes, elems, bottom, right, top, inner)

    print("Nodes:", len(nodes))
    print("Elements:", len(elems))