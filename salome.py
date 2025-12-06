import salome
import GEOM
from salome.geom import geomBuilder
import numpy as np
import sys
import importlib

sys.path.append('D:\\Airships\\Salome')

# This is to reload the local modules once they are changed. 
import plotter
importlib.reload(plotter)

salome.salome_init()
geom = geomBuilder.New()

# Development input paramaters
envelope_shape = "LOTTE"        # Shape of envelope
envelope_length = 100           # Axial length of envelope
envelope_resolution = 100       # Number of points taken axially

lobe_number = 3                 # Number of lobes of envelope (< 3)
lobe_offset_y = 6               # Distance between each lobe
lobe_offset_x = 13.333          # Central lobe offset (for 3 lobe design)
lobe_offset_z = 6               # Central lobe offset (for 3 lobe design)
central_lobe_length = 80        # Length of central lobe (for 3 lobe design)
central_lobe_shape = "LOTTE"    # Central lobe shape having a different shape.

fin_axial_offset = 80           # Distance between the leading edge of root chord of the fin and the nose. (< length - rc)
fin_thickness = 12              # Thickness-to-chord ratio * 100
fin_rc_length = 8               # Root chord length
fin_section_resolution = 40     # Number of points taken axially for the airfoil cross section
fin_taper_ratio = 0.5           # Ratio of tip chord length to root chord length
fin_height = 5                  # Perpendicular distance between root chord and the tip chord
fin_sweep_angle = 0             # Angle of sweep with respect to the normal to the cross section.
fin_number = 4                  # Number of fins. Should be even if there are more than one lobes.
fin_theta_pos = [0, 90]         # Angle positions for the fin with respect to vertical axis at top. In case of more than one lobe, this will only represent the theta positions of only one side of half of the fins.

# Modelling of envelope 

x_axis = geom.MakeVectorDXDYDZ(1, 0, 0)

def translate_object (object, x_offset, y_offset, z_offset):
    return geom.MakeTranslationTwoPoints(object, geom.MakeVertex(0, 0, 0), geom.MakeVertex(x_offset * lobe_offset_x, y_offset * lobe_offset_y, z_offset * lobe_offset_z))

def create_envelope (shape, length):
    gertler = plotter.GertlerEnvelope(plotter.STANDARD_ENVELOPES[shape], envelope_resolution, length)
    envelope_vertices = [geom.MakeVertex(x, y, 0) for x, y in gertler.points()]
    envelope_face = geom.MakeFace(geom.MakeWire([geom.MakeInterpol(envelope_vertices, False, False), geom.MakeLineTwoPnt(geom.MakeVertex(length, 0, 0), geom.MakeVertex(0, 0, 0))], 1e-7), 1)
    envelope = geom.MakeRevolution(envelope_face, x_axis, 2 * np.pi)
    
    return gertler, envelope

main_gertler, envelope = create_envelope(envelope_shape, envelope_length)
envelopes = [envelope] if lobe_number == 1 else [translate_object(envelope, 0, -1, 0), translate_object(envelope, 0, 1, 0)]

if lobe_number == 3:
    central_lobe = envelope if central_lobe_length == envelope_length else create_envelope(central_lobe_shape, central_lobe_length)[1]
    envelopes.append(translate_object(central_lobe, 1, 0, 1))

# Modelling of Fins

fin_axial_offset = min(fin_axial_offset, envelope_length - fin_rc_length)
fin_radial_offset = main_gertler.at(fin_axial_offset)
rc_axial_offset = fin_axial_offset
tc_axial_offset = fin_rc_length/2 * (1 - fin_taper_ratio) + fin_height * np.tan(np.radians(fin_sweep_angle))
rc_vertices = []
tc_vertices = []

for x, y in plotter.naca_airfoil_points(fin_thickness, fin_section_resolution, fin_rc_length):
    rc_vertices.append(geom.MakeVertex(fin_axial_offset + x, y, fin_radial_offset))
    tc_vertices.append(geom.MakeVertex(fin_axial_offset + tc_axial_offset + x * fin_taper_ratio, y * fin_taper_ratio, fin_radial_offset + fin_height))

rc_wire = geom.MakePolyline(rc_vertices, True)
tc_wire = geom.MakePolyline(tc_vertices, True)

rc_face = geom.MakeFace(rc_wire, True)
tc_face = geom.MakeFace(tc_wire, True)

midchord_direction = [geom.MakeVertex(0, 0, 0), geom.MakeVertex(tc_axial_offset, 0, fin_height)]
planform_surface = geom.MakePipeWithDifferentSectionsBySteps(
    [rc_wire, tc_wire],
    midchord_direction, 
    geom.MakePolyline(midchord_direction, False), 
)

x2, z2 = main_gertler.get_trailing_edge_intercept(fin_axial_offset, fin_rc_length)
fin = geom.MakeSolid(geom.MakeShell([planform_surface, rc_face, tc_face]))
fin = geom.MakeRotationThreePoints(
    fin, 
    geom.MakeVertex(fin_axial_offset, 0, fin_radial_offset),
    geom.MakeVertex(fin_axial_offset + fin_rc_length, 0, fin_radial_offset),
    geom.MakeVertex(x2, 0, z2)
)

fins = []

if lobe_number == 1:
    if not fin_theta_pos:
        fin_theta_pos = [i * (360 / fin_number) for i in range(0, fin_number)]

    fins = [geom.MakeRotation(fin, x_axis, -np.pi/2 + np.radians(theta)) for theta in fin_theta_pos]
else:
    for theta in fin_theta_pos:
        fins.append(translate_object(geom.MakeRotation(fin, x_axis, np.radians(theta)), 0, -1, 0))
        fins.append(translate_object(geom.MakeRotation(fin, x_axis, np.radians(-theta)), 0, 1, 0))

airship = geom.MakeFuseList(envelopes + fins)
airship_id = geom.addToStudy(airship, "Airship")

gg = salome.ImportComponentGUI("GEOM")
gg.createAndDisplayGO(airship_id)
gg.setDisplayMode(airship_id, 1)

# geom.ExportBREP(airship, 'D:\\Airships\\Salome\\test2.brep')

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()