# All the input paramaters required to generate the airship model.
# envelope_coeffs = (1.1518, -5.69072152, 27.47050632, -61.83093344, 58.45423403, -19.55488538)
envelope_length = 100           
envelope_coeffs = (1.2, -0.8779, -3.1206, 5.9936, -4.9138, 1.7187)
# envelope_diameter = 25.62853
envelope_diameter = 25.91344908
envelope_resolution = 100

lobe_number = 3
lobe_offset_x = 13.333
lobe_offset_y = 13.333/2
lobe_offset_z = 7        

central_lobe_coeffs = (1.1518, -5.69072152, 27.47050632, -61.83093344, 58.45423403, -19.55488538)
central_lobe_length = 80      
central_lobe_diameter = 20.5028319

fin_axial_offset = 80         
fin_thickness = 12            
fin_rc_length = 8             
fin_section_resolution = 40   
fin_taper_ratio = 0.5         
fin_height = 5                
fin_sweep_angle = 0           
fin_tip_angle = 0             
fin_number = 4                
fin_theta_pos = [0, 90]    

sheet_length_ratio = 0.75

directory_path = "D:\\Airships\\Salome"
final_object_name = "Airship"

# Import all the required modules
import salome
import GEOM
from salome.geom import geomBuilder
import numpy as np
import sys
import importlib

# Salome executes the python script files from its own directory so to import local modules, we have to 
# add our own path manually.
sys.path.append(directory_path)

# This is to reload the local modules once they are changed. For some reason, Salome does not reload the
# modules automatically.
# 
# NOTE: This is only a development thing and is to be removed when used in production.
import plotter
importlib.reload(plotter)

# Inititating the Geometry Module of Salome.
salome.salome_init()
geompy = geomBuilder.New()

# ---
# Elementary objects, directions and functions.
# ---

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0) 
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

# To translate any object to a position with offsets in multiples of lobe offsets.
def translate_object (object, x_offset, y_offset, z_offset):
    return geompy.MakeTranslationTwoPoints(object, O, geompy.MakeVertex(x_offset * lobe_offset_x, y_offset * lobe_offset_y, z_offset * lobe_offset_z))

# ---
# Modelling of envelope
# ---

# A function to create the basic envelope geometry along OX.
def create_envelope (coeffs, length, diameter):
    gertler = plotter.GertlerEnvelope(coeffs, length, diameter, envelope_resolution)
    envelope_vertices = [geompy.MakeVertex(x, y, 0) for x, y in gertler.points()]

    # NOTE: Do not make use of Polyline anywhere else as it results in revolution of solids
    # which do not undergo boolean operation properly. Make use of interpol and wires instead.
    # 
    # Given that airship hulls are hollow with some thickness, the revolution here probably
    # produces a filled solid (revolving a face). If so, we can later just extract the surface
    # of the formed solid for structural analysis. This is not an issue for CFD meshing atleast.
    #
    # Using interpol also results in a rough surface of revolution which is not the case when
    # using polyline. But, polyline fails in boolean operations. So, this is a tradeoff.
    #
    # TODO: Find a way to smoothen out the surface without any issues with boolean operations.
    envelope_wire = geompy.MakeWire([geompy.MakeInterpol(envelope_vertices, False, False), geompy.MakeLineTwoPnt(geompy.MakeVertex(length, 0, 0), O)], 1e-7)
    envelope_face = geompy.MakeFace(envelope_wire, 1)
    envelope = geompy.MakeRevolution(envelope_face, OX, 2 * np.pi)
    
    return gertler, envelope

# Creating the extreme lobe geometry first.
extreme_gertler, extreme_lobe = create_envelope(envelope_coeffs, envelope_length, envelope_diameter)

# If it is a single lobe design, there is only one extreme lobe which is the central lobe.
# If anything else (bi lobe or tri lobe), there are two extreme lobes whose translations are applied here.
lobes = [extreme_lobe] if lobe_number == 1 else [translate_object(extreme_lobe, 0, -1, 0), translate_object(extreme_lobe, 0, 1, 0)]

# If it is a tri lobe design, we have to model the central lobe as well.
if lobe_number == 3:
    # In case of no explicit specifications of central lobe coefficients, it means, the central
    # lobe is similar to the extreme lobes. If not, a new geometry has to be created.
    central_lobe = extreme_lobe if not central_lobe_coeffs else create_envelope(central_lobe_coeffs, central_lobe_length, central_lobe_diameter)[1]
    
    # Apply the necessary translations to the central lobe.
    lobes.append(translate_object(central_lobe, 1, 0, 1))

# ---
# Modelling of Fins
# ---

# Distance of the leading edge of root chord and tip chord from the axis of the envelope.
rc_radial_offset = extreme_gertler.at(fin_axial_offset)
tc_radial_offset = rc_radial_offset + fin_height

# Distance of the leading edge of root chord and tip chord from the nose of the envelope.
rc_axial_offset = fin_axial_offset
tc_axial_offset = rc_axial_offset + fin_rc_length/2 * (1 - fin_taper_ratio) + fin_height * np.tan(np.radians(fin_sweep_angle))

rc_vertices = []
tc_vertices = []

for x, y in plotter.naca_airfoil_points(fin_thickness, fin_section_resolution, fin_rc_length):
    rc_vertices.append(geompy.MakeVertex(rc_axial_offset + x, y, rc_radial_offset))
    # For tip chord, the points are linearly scaled by the taper ratio.
    tc_vertices.append(geompy.MakeVertex(tc_axial_offset + x * fin_taper_ratio, y * fin_taper_ratio, tc_radial_offset))

rc_wire = geompy.MakePolyline(rc_vertices, True)
tc_wire = geompy.MakePolyline(tc_vertices, True)

rc_face = geompy.MakeFace(rc_wire, True)
tc_face = geompy.MakeFace(tc_wire, True)

# Modelling the planform surface of the fin as pipe from root chord to tip chord along the mid chord line.
midchord_direction = [geompy.MakeVertex(rc_axial_offset, 0, 0), geompy.MakeVertex(tc_axial_offset, 0, fin_height)]
planform_surface = geompy.MakePipeWithDifferentSectionsBySteps(
    [rc_wire, tc_wire],
    midchord_direction, 
    geompy.MakePolyline(midchord_direction, False), 
)

# Rotating the fin so that the trailing edge of the root chord intersects with the extreme lobe surface.
x2, z2 = extreme_gertler.get_trailing_edge_intercept(rc_axial_offset, fin_rc_length)
fin = geompy.MakeSolid(geompy.MakeShell([planform_surface, rc_face, tc_face]))
fin = geompy.MakeRotationThreePoints(
    fin, 
    geompy.MakeVertex(rc_axial_offset, 0, rc_radial_offset),
    geompy.MakeVertex(rc_axial_offset + fin_rc_length, 0, rc_radial_offset),
    geompy.MakeVertex(x2, 0, z2)
)

fins = []

if lobe_number == 1:
    # If customised theta positions are not provided, they are evenly distributed.
    if not fin_theta_pos:
        fin_theta_pos = [i * (360 / fin_number) for i in range(0, fin_number)]

    # In case of single lobe, all fins lie on the same lobe.
    fins = [geompy.MakeRotation(fin, OX, -np.pi/2 + np.radians(theta)) for theta in fin_theta_pos]
else:
    # In case of bi lobe or tri lobe designs, fins are mirrored on both sides of the lobes.
    for theta in fin_theta_pos:
        fins.append(translate_object(geompy.MakeRotation(fin, OX, np.radians(theta)), 0, -1, 0))
        fins.append(translate_object(geompy.MakeRotation(fin, OX, np.radians(-theta)), 0, 1, 0))

# ---
# Modelling of Thin Fairings
# ---

fairings = []

# Creating a quadrilateral surface from the 4 points and give it a thickness in both directions normal to the surface.
def create_fairing_quad (p1, p2, p3, p4):
    fairing = geompy.MakeQuad4Vertices(geompy.MakeVertex(*p1), geompy.MakeVertex(*p2), geompy.MakeVertex(*p3), geompy.MakeVertex(*p4))
    normal = geompy.MakeVectorDXDYDZ(*np.cross(np.array(p2) - np.array(p1), np.array(p3) - np.array(p1)))
    fairings.append(geompy.MakePrismVecH2Ways(fairing, normal, 1e-7))

if sheet_length_ratio:
    sheet_length = envelope_length * sheet_length_ratio

    # In case of a bi lobe design, only one fairing sheet is required.
    if lobe_number == 2:
        create_fairing_quad(
            (envelope_length, -lobe_offset_y, 0),
            (envelope_length, lobe_offset_y, 0),
            (envelope_length - sheet_length, lobe_offset_y, 0),
            (envelope_length - sheet_length, -lobe_offset_y, 0)
        )

    # In case of a tri lobe design, two fairing sheets are required.
    elif lobe_number == 3:
        create_fairing_quad(
            (envelope_length, -lobe_offset_y, 0),
            (central_lobe_length + lobe_offset_x, 0, lobe_offset_z),
            (envelope_length - sheet_length, -lobe_offset_y, 0),
            (central_lobe_length + lobe_offset_x - sheet_length, 0, lobe_offset_z),
        )

        create_fairing_quad(
            (envelope_length, lobe_offset_y, 0),
            (central_lobe_length + lobe_offset_x, 0, lobe_offset_z),
            (envelope_length - sheet_length, lobe_offset_y, 0),
            (central_lobe_length + lobe_offset_x - sheet_length, 0, lobe_offset_z),
        )

# For the final airship model, fuse whichever compounds can be fused while making compound for the rest.
airship = geompy.MakeFuseList(lobes + fins)
airship = geompy.MakeCompound([airship] + fairings)             # Fusing fairings results in corrupted shapes.
airship_id = geompy.addToStudy(airship, final_object_name)

# If Salome GUI is present, display the final airship model and update the object browser.
if salome.sg.hasDesktop():
    gg = salome.ImportComponentGUI("GEOM")
    gg.createAndDisplayGO(airship_id)
    gg.setDisplayMode(airship_id, 1)

    salome.sg.updateObjBrowser()