import os
import subprocess
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
        self.lobe_number = self.paramaters["LOBE_NUMBER"]

        if self.lobe_number == 3:
            self.central_lobe = self.init_lobe("CENTRAL_LOBE", self.envelope)
        
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
    
    def run (self, open_gui = False, remove_temp_script = False):
        exec_script_path = os.path.join(DIR_PATH, EXEC_SCRIPT_FILE)
        open(exec_script_path, 'w').write(self.script())
        
        subprocess.run([self.salome_exec_path, "-g" if open_gui else "-t", exec_script_path], check=True)
        
        if remove_temp_script:
            os.remove(exec_script_path)