# All the input paramaters required to generate the airship model.
ENVELOPE_LENGTH = 100           
ENVELOPE_COEFFS = (1.2, -0.8779, -3.1206, 5.9936, -4.9138, 1.7187)
ENVELOPE_DIAMETER = 25.91344908
ENVELOPE_RESOLUTION = 100
ENVELOPE_TRUNCATION_RATIO = 0

LOBE_NUMBER = 3
LOBE_OFFSET_X = 13.333
LOBE_OFFSET_Y = 13.333/2
LOBE_OFFSET_Z = 7        

CENTRAL_LOBE_COEFFS = (1.1518, -5.69072152, 27.47050632, -61.83093344, 58.45423403, -19.55488538)
CENTRAL_LOBE_LENGTH = 80      
CENTRAL_LOBE_DIAMETER = 20.5028319

FIN_AXIAL_OFFSET = 80         
FIN_THICKNESS = 12            
FIN_RC_LENGTH = 8             
FIN_SECTION_RESOLUTION = 40   
FIN_TAPER_RATIO = 0.5         
FIN_HEIGHT = 5                
FIN_SWEEP_ANGLE = 0           
FIN_TIP_ANGLE = 10       
FIN_NUMBER = 4                
FIN_THETA_POS = [0, 90]    

SHEET_LENGTH_RATIO = 0.75

DIRECTORY_PATH = "D:\\Airships\\Salome"
FINAL_OBJECT_NAME = "Airship"

# Import all the required modules
import salome
import GEOM
from salome.geom import geomBuilder
import numpy as np
import sys
import importlib

# Salome executes the python script files from its own directory so to import local modules, we have to 
# add our own path manually.
sys.path.append(DIRECTORY_PATH)

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
    return geompy.MakeTranslationTwoPoints(object, O, geompy.MakeVertex(x_offset * LOBE_OFFSET_X, y_offset * LOBE_OFFSET_Y, z_offset * LOBE_OFFSET_Z))

# ---
# Modelling of envelope
# ---

# A function to create the basic envelope geometry along OX.
def create_envelope (coeffs, length, diameter):
    gertler = plotter.GertlerEnvelope(coeffs, length, diameter, ENVELOPE_RESOLUTION)
    envelope_vertices = [geompy.MakeVertex(x, y, 0) for x, y in gertler.points(ENVELOPE_TRUNCATION_RATIO)]

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
    envelope_edges = [geompy.MakeInterpol(envelope_vertices, False, False), geompy.MakeLineTwoPnt(geompy.MakeVertex(length * (1 - ENVELOPE_TRUNCATION_RATIO), 0, 0), O)]

    # When there is a truncation, we need to close the back end of the envelope as well. 
    if ENVELOPE_TRUNCATION_RATIO:
        # In case of truncation ratio being very small, the lost point of the envelope can end up being very
        # close to the axis resulting in errors.
        try: envelope_edges.append(geompy.MakeLineTwoPnt(geompy.MakeVertex(length * (1 - ENVELOPE_TRUNCATION_RATIO), geompy.PointCoordinates(envelope_vertices[-1])[1], 0), geompy.MakeVertex(length * (1 - ENVELOPE_TRUNCATION_RATIO), 0, 0)))
        except: pass

    envelope_wire = geompy.MakeWire(envelope_edges, 1e-7)
    envelope_face = geompy.MakeFace(envelope_wire, 1)
    envelope = geompy.MakeRevolution(envelope_face, OX, 2 * np.pi)
    
    return gertler, envelope

# Creating the extreme lobe geometry first.
extreme_gertler, extreme_lobe = create_envelope(ENVELOPE_COEFFS, ENVELOPE_LENGTH, ENVELOPE_DIAMETER)

# If it is a single lobe design, there is only one extreme lobe which is the central lobe.
# If anything else (bi lobe or tri lobe), there are two extreme lobes whose translations are applied here.
lobes = [extreme_lobe] if LOBE_NUMBER == 1 else [translate_object(extreme_lobe, 0, -1, 0), translate_object(extreme_lobe, 0, 1, 0)]

# If it is a tri lobe design, we have to model the central lobe as well.
if LOBE_NUMBER == 3:
    # In case of no explicit specifications of central lobe coefficients, it means, the central
    # lobe is similar to the extreme lobes. If not, a new geometry has to be created.
    central_lobe = extreme_lobe if not CENTRAL_LOBE_COEFFS else create_envelope(CENTRAL_LOBE_COEFFS, CENTRAL_LOBE_LENGTH, CENTRAL_LOBE_DIAMETER)[1]
    
    # Apply the necessary translations to the central lobe.
    lobes.append(translate_object(central_lobe, 1, 0, 1))

# ---
# Modelling of Fins
# ---

# Distance of the leading edge of root chord and tip chord from the axis of the envelope.
RC_RADIAL_OFFSET = extreme_gertler.at(FIN_AXIAL_OFFSET)
TC_RADIAL_OFFSET = RC_RADIAL_OFFSET + FIN_HEIGHT

# Distance of the leading edge of root chord and tip chord from the nose of the envelope.
RC_AXIAL_OFFSET = FIN_AXIAL_OFFSET
TC_AXIAL_OFFSET = RC_AXIAL_OFFSET + FIN_RC_LENGTH/2 * (1 - FIN_TAPER_RATIO) + FIN_HEIGHT * np.tan(np.radians(FIN_SWEEP_ANGLE))

# Scaling factors for fin tip chord due to tip angle rotation
COS_TIP_ANGLE = np.cos(np.radians(FIN_TIP_ANGLE))
SIN_TIP_ANGLE = np.sin(np.radians(FIN_TIP_ANGLE))

rc_vertices = []
tc_vertices = []

for x, y in plotter.naca_airfoil_points(FIN_THICKNESS, FIN_SECTION_RESOLUTION, FIN_RC_LENGTH):
    rc_vertices.append(geompy.MakeVertex(RC_AXIAL_OFFSET + x, y, RC_RADIAL_OFFSET))
    # For tip chord, the points are linearly scaled by the taper ratio.
    tc_vertices.append(geompy.MakeVertex(TC_AXIAL_OFFSET + x * FIN_TAPER_RATIO * COS_TIP_ANGLE, y * FIN_TAPER_RATIO, TC_RADIAL_OFFSET - x * FIN_TAPER_RATIO * SIN_TIP_ANGLE))

rc_wire = geompy.MakePolyline(rc_vertices, True)
tc_wire = geompy.MakePolyline(tc_vertices, True)

rc_face = geompy.MakeFace(rc_wire, True)
tc_face = geompy.MakeFace(tc_wire, True)

# Modelling the planform surface of the fin as pipe from root chord to tip chord along the mid chord line.
midchord_direction = [geompy.MakeVertex(RC_AXIAL_OFFSET, 0, 0), geompy.MakeVertex(TC_AXIAL_OFFSET, 0, FIN_HEIGHT)]
planform_surface = geompy.MakePipeWithDifferentSectionsBySteps(
    [rc_wire, tc_wire],
    midchord_direction, 
    geompy.MakePolyline(midchord_direction, False), 
)

# Rotating the fin so that the trailing edge of the root chord intersects with the extreme lobe surface.
x2, z2 = extreme_gertler.get_trailing_edge_intercept(RC_AXIAL_OFFSET, FIN_RC_LENGTH)
fin = geompy.MakeSolid(geompy.MakeShell([planform_surface, rc_face, tc_face]))
fin = geompy.MakeRotationThreePoints(
    fin, 
    geompy.MakeVertex(RC_AXIAL_OFFSET, 0, RC_RADIAL_OFFSET),
    geompy.MakeVertex(RC_AXIAL_OFFSET + FIN_RC_LENGTH, 0, RC_RADIAL_OFFSET),
    geompy.MakeVertex(x2, 0, z2)
)

fins = []

if LOBE_NUMBER == 1:
    # If customised theta positions are not provided, they are evenly distributed.
    if not FIN_THETA_POS:
        FIN_THETA_POS = [i * (360 / FIN_NUMBER) for i in range(0, FIN_NUMBER)]

    # In case of single lobe, all fins lie on the same lobe.
    fins = [geompy.MakeRotation(fin, OX, -np.pi/2 + np.radians(theta)) for theta in FIN_THETA_POS]
else:
    # In case of bi lobe or tri lobe designs, fins are mirrored on both sides of the lobes.
    for theta in FIN_THETA_POS:
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

if SHEET_LENGTH_RATIO:
    SHEET_LENGTH = ENVELOPE_LENGTH * SHEET_LENGTH_RATIO

    # In case of a bi lobe design, only one fairing sheet is required.
    if LOBE_NUMBER == 2:
        create_fairing_quad(
            (ENVELOPE_LENGTH, -LOBE_OFFSET_Y, 0),
            (ENVELOPE_LENGTH, LOBE_OFFSET_Y, 0),
            (ENVELOPE_LENGTH - SHEET_LENGTH, LOBE_OFFSET_Y, 0),
            (ENVELOPE_LENGTH - SHEET_LENGTH, -LOBE_OFFSET_Y, 0)
        )

    # In case of a tri lobe design, two fairing sheets are required.
    elif LOBE_NUMBER == 3:
        create_fairing_quad(
            (ENVELOPE_LENGTH, -LOBE_OFFSET_Y, 0),
            (CENTRAL_LOBE_LENGTH + LOBE_OFFSET_X, 0, LOBE_OFFSET_Z),
            (ENVELOPE_LENGTH - SHEET_LENGTH, -LOBE_OFFSET_Y, 0),
            (CENTRAL_LOBE_LENGTH + LOBE_OFFSET_X - SHEET_LENGTH, 0, LOBE_OFFSET_Z),
        )

        create_fairing_quad(
            (ENVELOPE_LENGTH, LOBE_OFFSET_Y, 0),
            (CENTRAL_LOBE_LENGTH + LOBE_OFFSET_X, 0, LOBE_OFFSET_Z),
            (ENVELOPE_LENGTH - SHEET_LENGTH, LOBE_OFFSET_Y, 0),
            (CENTRAL_LOBE_LENGTH + LOBE_OFFSET_X - SHEET_LENGTH, 0, LOBE_OFFSET_Z),
        )

# For the final airship model, fuse whichever compounds can be fused while making compound for the rest.
airship = geompy.MakeFuseList(lobes + fins)
airship = geompy.MakeCompound([airship] + fairings)             # Fusing fairings results in corrupted shapes.
airship_id = geompy.addToStudy(airship, FINAL_OBJECT_NAME)

# If Salome GUI is present, display the final airship model and update the object browser.
if salome.sg.hasDesktop():
    gg = salome.ImportComponentGUI("GEOM")
    gg.createAndDisplayGO(airship_id)
    gg.setDisplayMode(airship_id, 1)

    salome.sg.updateObjBrowser()