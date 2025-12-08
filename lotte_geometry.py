# geometry_handler.py (Modified Version: Full Intact Hull, Exact Width Sheet, NO FINS)

import numpy as np
import os
import tempfile
import math
import subprocess
from typing import Dict, Any

# --- GLOBAL CONSTANTS ---
SALOME_EXECUTABLE_PATH = "C:\\SALOME-9.15.0\\run_SALOME.bat"
NUM_POINTS_X_HIGH_RES = 1000
NUM_POINTS_C_DEFAULT = 11
NUM_ELLIPSE_POINTS = 50

# --- GEOMETRY PARAMETERS ---
LOBE_DIAMETER_RATIO = 0.35      # Lobe separation setting
STERN_BLEND_START_RATIO = 0.90
LOBE_TRUNCATION_RAT = 0.999     # Hull generated almost full length (No truncation)
SHEET_START_RATIO = 0.75        # Where the flat sheet STARTS (75% length)
TAIL_END_RATIO = 0.99
TAIL_TIP_WIDTH_MM = 5.0
TAIL_TIP_WIDTH_M = TAIL_TIP_WIDTH_MM / 1000.0

# --- FIN PARAMETERS (ALL REMOVED/DELETED) ---
# FIN_LENGTH_RATIO = 0.20
# FIN_HEIGHT_RATIO = 0.15
# FIN_THICKNESS_RATIO = 0.01
# FIN_OFFSET_RATIO = 0.10
# FIN_ANGLE_DEG_START = 45.0

# ---------------------------------------------------------------------------------
# --- PROFILE CALCULATION FUNCTIONS -----------------------------------------------
# ---------------------------------------------------------------------------------
def calculate_gertler_profile(params, L, x_norm):
    m_val, r0_val, r1_val, cp_val, l2d_val = params["m"], params["r0"], params["r1"], params["cp"], params["l2d"]
    A = np.array([[1, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1],
                  [m_val, m_val**2, m_val**3, m_val**4, m_val**5, m_val**6],
                  [1, 2*m_val, 3*m_val**2, 4*m_val**3, 5*m_val**4, 6*m_val**5],
                  [1, 2, 3, 4, 5, 6], [1/2, 1/3, 1/4, 1/5, 1/6, 1/7]])
    B = np.array([2*r0_val, 0, 1/4, 0, -2*r1_val, 1/4*cp_val])
    try:
        G = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        return None, None
    D = L / l2d_val
    y_sq_term = (G[0]*x_norm + G[1]*x_norm**2 + G[2]*x_norm**3 + G[3]*x_norm**4 + G[4]*x_norm**5 + G[5]*x_norm**6)
    y_sq_term[y_sq_term < 0] = 0
    R_x = D * np.sqrt(y_sq_term)
    return R_x, D

# ---------------------------------------------------------------------------------
# --- SALOME SCRIPT GENERATION ----------------------------------------------------
# ---------------------------------------------------------------------------------
def generate_salome_script(shape_name, params, hull_length):
    geom_type = params["type"]
    m_val, r0_val, r1_val, cp_val, l2d_val = params["m"], params["r0"], params["r1"], params["cp"], params["l2d"]

    geom_calc_func_base = f"""
import numpy as np 
# Gertler Profile Calculation Parameters
m, r0, r1, cp, l2d = {m_val}, {r0_val}, {r1_val}, {cp_val}, {l2d_val}
AIRSHIP_LENGTH = {hull_length}
LOBE_DIAMETER_RATIO = {LOBE_DIAMETER_RATIO}
STERN_BLEND_START_RATIO = {STERN_BLEND_START_RATIO} 
LOBE_TRUNCATION_RAT = {LOBE_TRUNCATION_RAT} 
SHEET_START_RATIO = {SHEET_START_RATIO}      
TAIL_END_RATIO = {TAIL_END_RATIO}
TAIL_TIP_WIDTH_M = {TAIL_TIP_WIDTH_M} 
NUM_ELLIPSE_POINTS = {NUM_ELLIPSE_POINTS}

AIRSHIP_DIAMETER = AIRSHIP_LENGTH / l2d
AIRSHIP_RADIUS = AIRSHIP_DIAMETER / 2.0
LOBE_HALF_OFFSET = (AIRSHIP_DIAMETER * LOBE_DIAMETER_RATIO) / 2.0 

# --- Profile Matrix Solution ---
A = np.array([
    [1, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1],
    [m, m**2, m**3, m**4, m**5, m**6],
    [1, 2*m, 3*m**2, 4*m**3, 5*m**4, 6*m**5],
    [1, 2, 3, 4, 5, 6], [1/2, 1/3, 1/4, 1/5, 1/6, 1/7]
])
B = np.array([2*r0, 0, 1/4, 0, -2*r1, 1/4*cp])
G = np.linalg.solve(A, B)

# --- Profile Calculation ---
X_profile_full = np.linspace(0, AIRSHIP_LENGTH, {NUM_POINTS_X_HIGH_RES})
x_norm_full = X_profile_full / AIRSHIP_LENGTH
y_sq_term = (G[0]*x_norm_full + G[1]*x_norm_full**2 + G[2]*x_norm_full**3 + G[3]*x_norm_full**4 + G[4]*x_norm_full**5 + G[5]*x_norm_full**6)
y_sq_term[y_sq_term < 0] = 0
R_profile_full = AIRSHIP_DIAMETER * np.sqrt(y_sq_term)

# --- Hull Truncation (Using LOBE_TRUNCATION_RAT=0.999 to keep the hull almost full length) ---
X_TRUNCATION_FINAL = AIRSHIP_LENGTH * LOBE_TRUNCATION_RAT
index_trunc_final = np.argmin(np.abs(X_profile_full - X_TRUNCATION_FINAL))
X_profile_hull = X_profile_full[:index_trunc_final]
R_profile_hull = R_profile_full[:index_trunc_final]
R_TRUNCATED_FINAL = R_profile_hull[-1]
X_profile_hull = np.append(X_profile_hull, X_TRUNCATION_FINAL)
R_profile_hull = np.append(R_profile_hull, R_TRUNCATED_FINAL)
R_profile_hull[0] = 0.0

R_MAX = AIRSHIP_RADIUS 

# Helper function for robust solid creation (used for mono-hull base)
def create_mono_hull(X_profile, R_profile):
    points_top = []
    # Use the profile array
    for i in range(len(X_profile)):
        points_top.append(geompy.MakeVertex(X_profile[i], R_profile[i], 0.0))
    points_bottom = []
    for i in range(len(X_profile) - 2, 0, -1): 
        points_bottom.append(geompy.MakeVertex(X_profile[i], 0.0, 0.0))
    all_points = points_top + points_bottom
    profile_wire = geompy.MakePolyline(all_points, True) 
    profile_face = geompy.MakeFace(profile_wire, 1)
    hull_solid = geompy.MakeRevolution(profile_face, OX, 2 * math.pi)
    return hull_solid

"""

    # --- GEOMETRY CONSTRUCTION ---
    if geom_type == "gertler":
        geometry_construction = """
# --- Monocoque Hull Construction (Single Lobe) ---
hull_shape = create_mono_hull(X_profile_hull, R_profile_hull)
"""

    elif geom_type == "gertler_bilobe":
        geometry_construction = """
# --- 1. Bi-Lobe Hull Construction (Full, Intact Shape) ---

# Create the base mono-hull solid (using the nearly full profile)
base_lobe_solid = create_mono_hull(X_profile_hull, R_profile_hull)

# Translate and Fuse the two solid lobes
lobe_A = geompy.MakeTranslation(base_lobe_solid, 0,  LOBE_HALF_OFFSET, 0)
lobe_B = geompy.MakeTranslation(base_lobe_solid, 0, -LOBE_HALF_OFFSET, 0)

hull_shape = geompy.MakeFuseList([lobe_A, lobe_B])


# -----------------------------------------------------------------
# --- 2. CREATE THIN RECTANGULAR SOLID FAIRING (Intersecting Sheet) ---
# -----------------------------------------------------------------

X_START_SHEET = AIRSHIP_LENGTH * SHEET_START_RATIO # Starts at 75% length
X_END_SHEET = AIRSHIP_LENGTH

# Dimensions of the box
DX_SHEET = X_END_SHEET - X_START_SHEET # Length of the fairing
DY_SHEET = (LOBE_HALF_OFFSET * 2.0) * 1.001 # Exact width limited to lobe separation + safety margin
DZ_SHEET = 1.0 / 1000000.0             # Height: Minimal Z thickness

# Define the lower-left corner for the box (MIN_X, MIN_Y, MIN_Z)
MIN_X_SHEET = X_START_SHEET
MIN_Y_SHEET = -DY_SHEET / 2.0  
MIN_Z_SHEET = -DZ_SHEET / 2.0 

thin_box_fairing = geompy.MakeBoxDXDYDZ(DX_SHEET, DY_SHEET, DZ_SHEET)
fairing_solid = geompy.MakeTranslation(thin_box_fairing, MIN_X_SHEET, MIN_Y_SHEET, MIN_Z_SHEET)
    
if fairing_solid is None:
    raise Exception("Failed to create Thin Box Fairing.")


# -----------------------------------------------------------------
# --- 3. FIN CONSTRUCTION REMOVED ---
# -----------------------------------------------------------------
# The fin geometry generation block was removed here.

"""

    else:
        raise ValueError("Unknown geometry type.")

    # --- Full SALOME Script Assembly ---
    script_content = f"""
import salome
salome.salome_init()
from salome.geom import geomBuilder
import GEOM
import math
import numpy as np
import os, sys 

geompy = geomBuilder.New()

O  = geompy.MakeVertex(0, 0, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)

{geom_calc_func_base}

{geometry_construction}

if hull_shape is None:
    raise Exception("Hull shape generation failed.")

# Add all separate solids to the study
geompy.addToStudy(hull_shape, "Full_Intact_BiLobe_Hull")
geompy.addToStudy(fairing_solid, "Intersecting_Flat_Sheet_Exact_Width")
# Fin objects were previously added here, but are now removed.

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()
"""
    return script_content

# ---------------------------------------------------------------------------------
# --- COORDINATE CALCULATION (Unchanged) ------------------------------------------
# ---------------------------------------------------------------------------------
def calculate_petal_coordinates(params, hull_length, num_petals, num_points_x=NUM_POINTS_X_HIGH_RES, num_points_c=NUM_POINTS_C_DEFAULT):
    geom_type = params["type"]
    L = hull_length
    X = np.linspace(0, L, num_points_x)
    x_norm = X / L
    if geom_type.startswith("gertler"): R_x, D_val = calculate_gertler_profile(params, L, x_norm)
    else: return [], 0, 0
    if R_x is None: return [], 0, 0
    delta_phi = 2 * math.pi / num_petals
    phi_half_width = delta_phi / 2.0
    Phi = np.linspace(-phi_half_width, phi_half_width, num_points_c)
    coords_2D = []
    for i in range(num_points_x):
        r = R_x[i]; x = X[i]
        for j in range(num_points_c):
            phi = Phi[j]
            C = r * phi
            coords_2D.append((x, C))
    return coords_2D, num_points_x, num_points_c

# ---------------------------------------------------------------------------------
# --- MAIN EXECUTION FUNCTION -----------------------------------------------------
# ---------------------------------------------------------------------------------
def generate_and_run_geometry(shape_name, params, hull_length, num_petals, dat_filename):
    salome_status_msg = ""; dat_status_msg = "";

    # 1. Generate and run SALOME script
    try:
        if not os.path.exists(SALOME_EXECUTABLE_PATH):
            raise FileNotFoundError(f"SALOME launch script not found at:\n{SALOME_EXECUTABLE_PATH}")

        script_content = generate_salome_script(shape_name, params, hull_length)
        script_filename = os.path.join(tempfile.gettempdir(), "salome_airship_script.py")

        with open(script_filename, "w") as f:
            f.write(script_content)

        command = [SALOME_EXECUTABLE_PATH, script_filename]
        # Starting the SALOME process without waiting for it to finish
        subprocess.Popen(command)

        # There are now only two objects added to the study: the hull and the sheet.
        salome_status_msg = (f"SALOME launched with {shape_name} **(MODIFIED GEOMETRY)**. Check your taskbar. You will see two separate objects in the study.")

    except Exception as e:
        salome_status_msg = f"Failed to execute SALOME:\n{e}"

    # 2. Generate DAT file coordinates
    if dat_filename:
        try:
            coords_2D, num_x, num_c = calculate_petal_coordinates(params, hull_length, num_petals)

            with open(dat_filename, "w") as f:
                f.write(f"{shape_name} 2D Developed Petal (X vs Circumferential Distance C)\n")
                f.write(f"{num_x}\t{num_c}\n")
                for x, c in coords_2D: f.write(f"{x: .6f}\t{c: .6f}\n")

            dat_status_msg = ("2D Developed Petal file successfully saved to:\n" f"{dat_filename}\n(File contains X vs. C coordinates)")

        except Exception as e:
            dat_status_msg = ("Error writing 2D coordinate file or during calculation:\n" f"{e}")
    else:
        dat_status_msg = "DAT file generation canceled by user."

    final_message = (f"--- SALOME Status ---\n{salome_status_msg}\n\n" f"--- DAT File Status ---\n{dat_status_msg}")
    return final_message