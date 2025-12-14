import numpy as np
from geometry import AirshipGeometry, plot_petal_profile, STANDARD_ENVELOPES, GertlerEnvelope

parameters = {
    "ENVELOPE_PARAMS": STANDARD_ENVELOPES["Wang"],
    # "ENVELOPE_LENGTH": 100,
    # "ENVELOPE_RESOLUTION": 100,
    "VOLUME": 500,

    "LOBE_NUMBER": 3,
    "LOBE_OFFSET_X": 13.333,
    "LOBE_OFFSET_Y": 13.333 / 2,
    "LOBE_OFFSET_Z": 7,

    "CENTRAL_LOBE_LENGTH": 80,

    "FIN_AXIAL_OFFSET": 80,
    "FIN_THICKNESS": 12,
    "FIN_RC_LENGTH": 8,
    "FIN_SECTION_RESOLUTION": 40,
    "FIN_TAPER_RATIO": 0.5,
    "FIN_HEIGHT": 5,
    "FIN_NUMBER": 4,
}

geometry = AirshipGeometry(parameters, "C:\\SALOME-9.15.0\\run_salome.bat")

print(geometry.volume())

# geometry.run_salome(open_gui=False, remove_temp_script=True, export_format='BREP')