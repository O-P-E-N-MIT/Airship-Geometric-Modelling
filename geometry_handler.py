# geometry_handler.py (FINAL VERSION with Tapered End Fix and High Resolution)

import os
import tempfile
import math
import numpy as np
import subprocess
import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt

# --- GLOBAL CONSTANTS ---
# NOTE: Update this path to match your SALOME installation location!
SALOME_EXECUTABLE_PATH = "C:\\SALOME-9.14.0\\run_SALOME.bat"
# Increased resolution for accurate tip plotting
NUM_POINTS_X_HIGH_RES = 1000
NUM_POINTS_C_DEFAULT = 11

# --- SALOME SCRIPT GENERATION (3D Model) ---

def generate_salome_script(shape_name, params, hull_length):
    """Generates the Python script content for SALOME."""
    m_val = params["m"]
    r0_val = params["r0"]
    r1_val = params["r1"]
    cp_val = params["cp"]
    l2d_val = params["l2d"]

    # The SALOME script content remains identical for 3D visualization.
    script_content = f"""
import salome
salome.salome_init()

from salome.geom import geomBuilder
import GEOM
import math
import numpy as np

geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)

def generate_hull_shape(m, r0, r1, cp, l2d, L=100.0, noe={NUM_POINTS_X_HIGH_RES}):

    m_val = m

    A = np.array([
        [1, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1],
        [m_val, m_val**2, m_val**3, m_val**4, m_val**5, m_val**6],
        [1, 2*m_val, 3*m_val**2, 4*m_val**3, 5*m_val**4, 6*m_val**5],
        [1, 2, 3, 4, 5, 6],
        [1/2, 1/3, 1/4, 1/5, 1/6, 1/7]
    ])

    B = np.array([2*r0, 0, 1/4, 0, -2*r1, 1/4*cp])

    try:
        G = np.linalg.solve(A, B)
        G = np.round(G, 4)
    except np.linalg.LinAlgError:
        print("Matrix is singular or ill-conditioned. Cannot solve for coefficients.")
        return None

    D = L / l2d

    X = np.linspace(0, L, noe)
    x_norm = X / L

    y_sq_term = (
        G[0]*x_norm + G[1]*x_norm**2 + G[2]*x_norm**3 +
        G[3]*x_norm**4 + G[4]*x_norm**5 + G[5]*x_norm**6
    )

    y_sq_term[y_sq_term < 0] = 0
    y = D * np.sqrt(y_sq_term)

    points = []
    for i in range(noe):
        points.append(geompy.MakeVertex(X[i], y[i], 0.0))

    top_wire = geompy.MakePolyline(points, False)
    # Hull Revolution 
    hull_solid = geompy.MakeRevolution(top_wire, OX, 2 * math.pi)

    return hull_solid

m = {m_val}
r0 = {r0_val}
r1 = {r1_val}
cp = {cp_val}
l2d = {l2d_val}
AIRSHIP_LENGTH = {hull_length}
AIRSHIP_RADIUS = AIRSHIP_LENGTH / l2d

# Using high resolution for 3D model generation too
hull_shape = generate_hull_shape(m, r0, r1, cp, l2d, L=AIRSHIP_LENGTH, noe={NUM_POINTS_X_HIGH_RES})

if hull_shape is None:
    raise Exception("Hull shape generation failed due to singular matrix.")

geompy.addToStudy(hull_shape, "Airship_Hull_{shape_name}".replace(" ", "_").replace("(", "").replace(")", ""))

b = AIRSHIP_RADIUS
a = AIRSHIP_LENGTH / 2.0

FIN_CHORD = a * 0.50
FIN_SPAN_EXPOSED = b * 0.8
FIN_THICKNESS_RATIO = 0.16
FIN_POSITION_X = AIRSHIP_LENGTH * 0.95

FIN_EMBEDDED_DEPTH = b * 0.50 + 10.0
TOTAL_FIN_LENGTH = FIN_SPAN_EXPOSED + FIN_EMBEDDED_DEPTH
FIN_NUM_POINTS = 50

Fin_Center_Axis = geompy.MakeVertex(FIN_POSITION_X, 0.0, 0.0)
X_TRANSLATION = FIN_POSITION_X - FIN_CHORD / 2.0
Z_ROOT_TRANSLATION = -TOTAL_FIN_LENGTH / 2.0

points_2D = []
x_coords = [(FIN_CHORD / 2.0) * (1.0 - math.cos(i * math.pi / (2*FIN_NUM_POINTS - 1))) for i in range(FIN_NUM_POINTS)]
def get_y_t(x_norm, t_c_ratio, chord_length):
    return 5.0 * t_c_ratio * chord_length * (
            0.2969 * math.sqrt(x_norm) - 0.1260 * x_norm -
            0.3516 * x_norm**2 + 0.2843 * x_norm**3 - 0.1015 * x_norm**4
    )

for x in x_coords:
    x_norm = x / FIN_CHORD
    y_t = get_y_t(x_norm, FIN_THICKNESS_RATIO, FIN_CHORD)
    points_2D.append(geompy.MakeVertex(x, y_t, 0.0))

for i in range(len(x_coords) - 2, 0, -1):
    x = x_coords[i]
    x_norm = x / FIN_CHORD
    y_t = get_y_t(x_norm, FIN_THICKNESS_RATIO, FIN_CHORD)
    points_2D.append(geompy.MakeVertex(x, -y_t, 0.0))

naca_wire = geompy.MakePolyline(points_2D, True)
naca_face = geompy.MakeFace(naca_wire, 1)
P1 = O
P2 = geompy.MakeVertex(0, 0, TOTAL_FIN_LENGTH)
fin_solid_proto = geompy.MakePrism(naca_face, P1, P2)

fin_translated = geompy.MakeTranslation(
    fin_solid_proto,
    X_TRANSLATION,
    0.0,
    Z_ROOT_TRANSLATION
)

# Rotation sequence to create four fins 
Fin_1_Vert = geompy.MakeRotation(fin_translated, OX, -math.pi / 2.0, Fin_Center_Axis)
geompy.addToStudy(Fin_1_Vert, "Fin_1_Vert_PosZ")

Fin_2_Horiz_posY = geompy.MakeRotation(Fin_1_Vert, OX, -math.pi / 2.0, Fin_Center_Axis)
geompy.addToStudy(Fin_2_Horiz_posY, "Fin_2_Horiz_posY")

Fin_3_Vert_negZ = geompy.MakeRotation(Fin_1_Vert, OX, -math.pi, Fin_Center_Axis)
geompy.addToStudy(Fin_3_Vert_negZ, "Fin_3_Vert_negZ")

Fin_4_Horiz_negY = geompy.MakeRotation(Fin_1_Vert, OX, math.pi / 2.0, Fin_Center_Axis)
geompy.addToStudy(Fin_4_Horiz_negY, "Fin_4_Horiz_negY")

# --- CRITICAL STEP: Combining all parts into one compound ---
all_parts = [hull_shape, Fin_1_Vert, Fin_2_Horiz_posY, Fin_3_Vert_negZ, Fin_4_Horiz_negY]
final_geometry = geompy.MakeCompound(all_parts)
geompy.addToStudy(final_geometry, "Combined_Airship_Compound_for_Meshing_V_Final")

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()
"""
    return script_content

# --- 2D PETAL CALCULATION (FIXED IMPLEMENTATION) ---

def calculate_petal_coordinates(params, hull_length, num_petals, num_points_x=NUM_POINTS_X_HIGH_RES, num_points_c=NUM_POINTS_C_DEFAULT):
    """
    Calculates the X (Longitudinal) and C (Circumferential Arc Length) coordinates
    for the 2D Developed Petal shape, mathematically enforcing zero width at the ends
    and using high resolution.
    """
    m_val = params["m"]
    r0_val = params["r0"]
    r1_val = params["r1"]
    cp_val = params["cp"]
    l2d_val = params["l2d"]
    L = hull_length

    # Matrix solver logic to get radius R(x)
    A = np.array([
        [1, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1], [m_val, m_val**2, m_val**3, m_val**4, m_val**5, m_val**6],
        [1, 2*m_val, 3*m_val**2, 4*m_val**3, 5*m_val**4, 6*m_val**5], [1, 2, 3, 4, 5, 6],
        [1/2, 1/3, 1/4, 1/5, 1/6, 1/7]
    ])
    B = np.array([2*r0_val, 0, 1/4, 0, -2*r1_val, 1/4*cp_val])

    try:
        G = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        print("Matrix is singular. Cannot calculate coordinates.")
        return [], 0, 0

    D = L / l2d_val
    X = np.linspace(0, L, num_points_x)
    x_norm = X / L

    # Calculate R(x)^2
    y_sq_term = (G[0]*x_norm + G[1]*x_norm**2 + G[2]*x_norm**3 +
                 G[3]*x_norm**4 + G[4]*x_norm**5 + G[5]*x_norm**6)
    y_sq_term[y_sq_term < 0] = 0
    R_x = D * np.sqrt(y_sq_term)

    # --- MATHEMATICAL CORRECTION: ENFORCING BOUNDARY CONDITIONS ---
    # This is the core logical correction: R(0)=0 and R(L)=0
    R_x[0] = 0.0
    R_x[-1] = 0.0
    # -----------------------------------------------------------------

    # --- PETAL DEVELOPMENT LOGIC ---
    delta_phi = 2 * math.pi / num_petals
    phi_half_width = delta_phi / 2.0
    Phi = np.linspace(-phi_half_width, phi_half_width, num_points_c)

    coords_2D = []

    # Iterate over axial points
    for i in range(num_points_x):
        r = R_x[i] # Current (corrected) radius
        x = X[i]   # Current axial position

        # Iterate over circumferential points (C)
        for j in range(num_points_c):
            phi = Phi[j]
            # C = R(X) * phi
            C = r * phi

            coords_2D.append((x, C))

    return coords_2D, num_points_x, num_points_c

# --- FUNCTION FOR PLOTTING AND SAVING PNG ---

def plot_and_save_profile(coords_2D, shape_name, dat_filename, num_points_x, num_points_c, num_petals):
    """Plots the 2D Developed Petal (X vs Circumferential Distance C) and saves it as a PNG."""
    try:
        # Reshape the coordinates back into a grid for plotting the boundaries
        X_flat = [c[0] for c in coords_2D]
        C_flat = [c[1] for c in coords_2D]
        X_grid = np.array(X_flat).reshape((num_points_x, num_points_c))
        C_grid = np.array(C_flat).reshape((num_points_x, num_points_c))

        png_filename = os.path.splitext(dat_filename)[0] + ".png"

        plt.figure(figsize=(10, 8))

        # Plot the side edges (longitudinal lines) - these define the petal shape
        plt.plot(X_grid[:, 0], C_grid[:, 0], 'r-', label='Longitudinal Edge (Boundary)', linewidth=2)
        plt.plot(X_grid[:, -1], C_grid[:, -1], 'r-', linewidth=2)

        # Plot the centerline (where C=0)
        plt.plot(X_grid[:, num_points_c // 2], C_grid[:, num_points_c // 2], 'k--', label='Petal Centerline', linewidth=1)

        # Plot the axial lines (cross-sections) for start and end
        # Since R(0)=0 and R(L)=0, these lines now correctly collapse to points at the origin/terminus.
        plt.plot(X_grid[0, :], C_grid[0, :], 'b-', label='Axial Boundary', linewidth=2)
        plt.plot(X_grid[num_points_x-1, :], C_grid[num_points_x-1, :], 'b-', linewidth=2)

        plt.title(f'2D Developed Petal: {shape_name} ({num_petals} Petals)', fontsize=16)
        plt.xlabel('Axial Position (X) [units]', fontsize=12)
        plt.ylabel('Circumferential Distance (C) [units]', fontsize=12)
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.legend()
        plt.axis('equal') # Crucial for showing the correct tapered shape

        plt.savefig(png_filename, dpi=300, bbox_inches='tight')
        plt.close()

        return f"\nPlot successfully saved to:\n{png_filename}"

    except Exception as e:
        return f"\nError generating plot:\n{e}"

# --- MAIN EXECUTION FUNCTION ---

def generate_and_run_geometry(shape_name, params, hull_length, num_petals, dat_filename):
    """Handles the creation of the SALOME script, launches SALOME, writes the 2D DAT file, and plots it."""

    salome_status_msg = ""
    plot_status_msg = ""
    coords_2D = []
    num_x = 0
    num_c = 0

    # 1. SALOME Script Generation and Execution (3D model generation)
    try:
        if not os.path.exists(SALOME_EXECUTABLE_PATH):
            raise FileNotFoundError(f"SALOME launch script not found at:\n{SALOME_EXECUTABLE_PATH}")

        script_content = generate_salome_script(shape_name, params, hull_length)
        script_filename = os.path.join(tempfile.gettempdir(), "salome_airship_script.py")

        with open(script_filename, "w") as f:
            f.write(script_content)

        command = [SALOME_EXECUTABLE_PATH, script_filename]
        subprocess.Popen(command)
        salome_status_msg = f"SALOME launched with {shape_name} 3D geometry."

    except Exception as e:
        salome_status_msg = f"Failed to execute SALOME:\n{e}"

    # 2. Calculate and Write 2D DAT File
    dat_status_msg = ""

    if dat_filename:
        try:
            # Calculation returns the list of points and the grid dimensions
            coords_2D, num_x, num_c = calculate_petal_coordinates(params, hull_length, num_petals)

            with open(dat_filename, "w") as f:
                f.write(f"{shape_name} 2D Developed Petal (X vs Circumferential Distance C)\n")
                # Write the grid dimensions (Nx, Nc)
                f.write(f"{num_x}\t{num_c}\n")

                # Write X and C (Circumferential Distance) coordinates
                for x, c in coords_2D:
                    f.write(f"{x: .6f}\t{c: .6f}\n")

            dat_status_msg = f"2D Developed Petal file successfully saved to:\n{dat_filename}\n(File contains X vs. C coordinates)"

            # 3. Plot the Profile
            plot_status_msg = plot_and_save_profile(coords_2D, shape_name, dat_filename, num_x, num_c, num_petals)

        except Exception as e:
            dat_status_msg = f"Error writing 2D coordinate file or during calculation:\n{e}"
    else:
        dat_status_msg = "DAT file generation canceled by user."

    # 4. Return Combined Status Message
    final_message = f"--- SALOME Status ---\n{salome_status_msg}\n\n--- DAT File Status ---\n{dat_status_msg}{plot_status_msg}"
    return final_message