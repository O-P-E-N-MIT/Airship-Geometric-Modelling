import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np

from plotter import GertlerEnvelope, STANDARD_ENVELOPES

# Main directory path where the executing script file is located
DIR_PATH = os.path.dirname(os.path.abspath(__file__))

# Base salome script file
BASE_SCRIPT_FILE = "salome.py"

# The temporary script file to be executed in Salome
EXEC_SCRIPT_FILE = "exec_salome_temp.py"

# The line in the base script file which indicates the end of input parameters
BREAKER = "# INPUT PARAMETERS END"

# Extract the base script content excluding the input parameters
BASE_SCRIPT = open(os.path.join(DIR_PATH, BASE_SCRIPT_FILE), 'r').read().split(BREAKER, 1)[1]

class AirshipGeometry:

    # TODO: Add ways to make optional paramaters.
    def __init__ (self, parameters, salome_exec_path):
        self.parameters = parameters
        self.salome_exec_path = salome_exec_path

        self.envelope = self.init_lobe("ENVELOPE")
        self.lobe_number = self.parameters["LOBE_NUMBER"]

        if self.lobe_number == 3:
            self.central_lobe = self.init_lobe("CENTRAL_LOBE")
        
        if "FINAL_OBJECT_NAME" not in self.parameters:
            self.parameters["FINAL_OBJECT_NAME"] = "Airship"
        
        self.parameters["DIRECTORY_PATH"] = DIR_PATH

    def init_lobe (self, id):
        len = self.parameters[f"{id}_LENGTH"]
        res = self.parameters[f"{id}_RESOLUTION"]
        envelope = None

        if f"{id}_SHAPE" in self.parameters:
            shape_name = self.parameters[f"{id}_SHAPE"]

            if shape_name in STANDARD_ENVELOPES:
                env_params = STANDARD_ENVELOPES[shape_name]
                envelope = GertlerEnvelope.from_parameters(env_params, len, res)
            
        elif f"{id}_PARAMS" in self.parameters:
            env_params = self.parameters[f"{id}_PARAMS"]
            envelope = GertlerEnvelope.from_parameters(env_params, len, res)
        
        elif f"{id}_COEFFS" in self.parameters:
            env_coeffs = self.parameters[f"{id}_COEFFS"]
            envelope = GertlerEnvelope(env_coeffs, len, self.parameters[f"{id}_DIAMETER"], res)

        self.parameters[f"{id}_COEFFS"] = [float(c) for c in envelope.coeffs]
        self.parameters[f"{id}_DIAMETER"] = envelope.diameter

        return envelope
    
    def script (self):
        param_lines = []
        for key, value in self.parameters.items():
            if isinstance(value, str):
                param_lines.append(f'{key} = "{value}"')
            else:
                param_lines.append(f"{key} = {value}")
        
        param_script = "\n".join(param_lines)
        full_script = f"{param_script}\n\n{BASE_SCRIPT}"
        return full_script
    
    def run_salome (self, open_gui = False, remove_temp_script = False):
        exec_script_path = os.path.join(DIR_PATH, EXEC_SCRIPT_FILE)
        open(exec_script_path, 'w').write(self.script())
        
        subprocess.run([self.salome_exec_path, "-g" if open_gui else "-t", exec_script_path], check=True)
        
        if remove_temp_script:
            os.remove(exec_script_path)


def plot_petal_profile (envelope, num_petals, num_points_c, dat_filename, shape_name='Envelope'):
    coords_2D = envelope.petal_coordinates(num_petals, num_points_c)
    num_points_x = envelope.n

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