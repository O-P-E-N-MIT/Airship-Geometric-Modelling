from geometry import AirshipGeometry, plot_petal_profile

parameters = {
    "ENVELOPE_LENGTH": 100,
    "ENVELOPE_SHAPE": "Wang",
    "ENVELOPE_RESOLUTION": 100,
    "ENVELOPE_TRUNCATION_RATIO": 0,

    "LOBE_NUMBER": 3,
    "LOBE_OFFSET_X": 13.333,
    "LOBE_OFFSET_Y": 13.333 / 2,
    "LOBE_OFFSET_Z": 7,

    "CENTRAL_LOBE_SHAPE": "LOTTE",
    "CENTRAL_LOBE_LENGTH": 80,
    "CENTRAL_LOBE_RESOLUTION": 100,

    "FIN_AXIAL_OFFSET": 80,
    "FIN_THICKNESS": 12,
    "FIN_RC_LENGTH": 8,
    "FIN_SECTION_RESOLUTION": 40,
    "FIN_TAPER_RATIO": 0.5,
    "FIN_HEIGHT": 5,
    "FIN_SWEEP_ANGLE": 0,
    "FIN_TIP_ANGLE": 10,
    "FIN_NUMBER": 4,
    "FIN_THETA_POS": [0, 90],

    "SHEET_LENGTH_RATIO": 0.75
}

geometry = AirshipGeometry(parameters, "C:\\SALOME-9.15.0\\run_salome.bat")

print(geometry.run_salome(open_gui=True, remove_temp_script=True))
plot_and_save_profile(geometry.envelope, 3, 100, "envelope_profile.dat", shape_name="Envelope")