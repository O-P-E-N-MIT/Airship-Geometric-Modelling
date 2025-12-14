# Input Variables for Gertler Envelope Generation
#
# m1 = Point of maximum thickness
# r0 = Nose radius at the front of the airship
# r1 = Tail radius at the back of the airship
# cp = Prismatic coefficient
# l2d = Length to diameter ratio
#
# Input Variables for a symmetric NACA 4 digit Airfoil Generation
#
# t = Maximum thickness as percentage of chord

import math
import numpy as np

from shapely.geometry import Point
from shapely.ops import unary_union

from scipy.optimize import fsolve

STANDARD_ENVELOPES = {
    "Sphere":   (0.5,   0.5,   0.5,   0.667, 1.000),
    "GNVR":     (0.415, 0.600, 0.180, 0.615, 3.044),
    "ZHIYUAN-1":(0.419, 0.337, 0.251, 0.651, 3.266),
    "Wang":     (0.404, 0.600, 0.100, 0.610, 3.859),
    "NPL":      (0.432, 0.589, 0.425, 0.667, 4.000),
    "LOTTE":    (0.4502, 0.5759, 0.10, 0.5170, 3.9019)
}

class GertlerEnvelope:

    # Input variables for Gertler Envelope
    #
    # coeffs = Coefficients of the Gertler polynomial
    # length = Length of the envelope
    # diameter = Maximum diameter of the envelope
    # n = Number of points to be generated along the length
    #
    # NOTE: Coefficients for Gertler Enevelope is to be of the following form. 
    # a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4 + a5 * x^5 + a6 * x^6
    def __init__ (self, coeffs, length, diameter, n):
        self.coeffs = list(coeffs)
        self.length = length
        self.diameter = diameter
        self.n = n

    # Get the radius of the envelope at a particular axial position x
    def at (self, x):
        x = x / self.length
        y = self.diameter * (self.coeffs[0]*x + self.coeffs[1]*x**2 + self.coeffs[2]*x**3 + self.coeffs[3]*x**4 + self.coeffs[4]*x**5 + self.coeffs[5]*x**6)**0.5
        return y
    
    # Returns an array of points representing the Gertler envelope.
    #
    # NOTE: These points are not linearly scalable with length like NACA airfoils.
    # NOTE: I have made use of numpy array here but an iterator function would be memory efficient I guess?
    def points (self, truncation = 0, tuples = True):
        X = np.linspace(0, 1 - truncation, self.n)
        R = np.polyval(self.ordered_coeffs, X)
        R[R < 0] = 0

        # In case if there is no truncation, the final point must lie on the axis.
        if not truncation:
            R[-1] = 0 

        R = self.diameter * np.sqrt(R)
        X = self.length * X

        # TODO: Zipping elements maybe a bad idea.
        return zip(X, R) if tuples else X, R
    
    def petal_coordinates (self, petal_number, circumferential_divisons):
        delta_phi = 2 * np.pi / petal_number
        phi_half_width = delta_phi / 2.0
        Phi = np.linspace(-phi_half_width, phi_half_width, circumferential_divisons)

        coords_2D = []

        for x, r in self.points():
            for j in range(circumferential_divisons):
                phi = Phi[j]
                # C = R(X) * phi
                C = r * phi

                coords_2D.append((x, C))

        return coords_2D
    
    # Returns the volume of the envelope (mono lobe).
    def volume (self):
        # pi x D^2 x L x cp/4
        return np.pi * (self.diameter**2) * self.length * np.dot(self.coeffs, [1/2, 1/3, 1/4, 1/5, 1/6, 1/7])
    
    # Returns the volume of a bilobed envelope.
    #
    # f = Distance of the extreme lobe from X axis.
    def volume_bilobe (self, f):
        # If there is no distance between lobes, both the lobes coincide.
        if f == 0:
            return self.volume()

        X, R = self.points(tuples=False)    # Array of coordinates of the envelope
        A = np.zeros(self.n)                # Array of intersections of cross section

        # Calculation intersecting cross sectional area
        r = R[f < R]
        A[f < R] =  2 * (r**2 * np.acos(f/r) - f * np.sqrt(r**2 - f**2))

        # 2 * (Volume of the lobes) - Intersection of both the lobes
        return 2 * self.volume() - np.trapezoid(A, X)
    
    # Returns the volume of a trilobed envelope.
    #
    # e = Distance of central lobe from Y axis
    # f = Distance of extreme lobe from X axis
    # g = Distance of central lobe from Z axis
    #
    # NOTE: For trilobe volume calculation, shapely module is used. Try to find some other way to compute
    # it other way which is more accurate and faster. Shapely calculations are not accurate but the error is
    # often negligible.
    #
    # TODO: Test this properly as I only tried only with few specific cases. Later, optimise the code.
    def volume_trilobe (self, e, f, g, central_lobe = None):
        # If central lobe is not specified, itself is taken as central lobe.
        central_lobe = central_lobe or self

        # Determining how long should the discretized X axis should be taken for calculation.
        last = max(self.length, e + central_lobe.length)
        # NOTE: Number of elements taking the same as the main envelope is not a good idea.
        X = np.linspace(0, last, self.n)

        # Getting the indices of position of start of central lobe, end of central lobe, end of extreme lobe in the
        # discretized X axis. This may not be the best way to do this.
        extreme_lobe_end_index = np.argmin(np.abs(X - self.length))
        central_lobe_start_index = np.argmin(np.abs(X - e))
        central_lobe_end_index = np.argmin(np.abs(X - (e + central_lobe.length)))

        # Due to discretization of X axis, the extreme lobe length and distance of central lobe from Y axis will be
        # different than the ideal values specified. Higher the value of n, closer will be to their ideal values.
        # TODO: Again, this may not be the best way to do it.
        extreme_lobe_length = X[extreme_lobe_end_index]     # Ideal value: self.envelope
        central_lobe_start = X[central_lobe_start_index]    # Ideal value: e

        # Radius of extreme lobe along the X axis
        # NOTE: This maybe a crude way to do this wihout using .points()
        R = np.zeros(self.n)
        R[:extreme_lobe_end_index+1] = np.polyval(self.ordered_coeffs, X[:extreme_lobe_end_index+1] / extreme_lobe_length)
        R[extreme_lobe_end_index:] = 0
        R = self.diameter * np.sqrt(R)

        # Radius of central lobe along the X axis.
        R2 = np.zeros(self.n)
        R2[central_lobe_start_index:] = np.polyval(central_lobe.ordered_coeffs, (X[central_lobe_start_index:] - central_lobe_start) / central_lobe.length)
        R2[central_lobe_end_index:] = 0
        R2 = central_lobe.diameter * np.sqrt(R2)

        # Array of cross section area along the axis.
        A = np.zeros(self.n)
        
        for i in range(self.n):
            r = R[i]
            r2 = R2[i]

            # Union of all the circle areas
            A[i] = unary_union([Point(-f, 0).buffer(r), Point(f, 0).buffer(r), Point(0, g).buffer(r2)]).area

        # Integrating the area of cross section to get the volume.
        return np.trapezoid(A, X)
    
    # Returns the coordinates of points on the envelope which intercepts the trailing edge of a fin.
    def get_trailing_edge_intercept (self, x, rc):
        y = self.at(x)
        chord_length = 0
        h = self.length / self.n
        x1 = x
        y1 = y

        while chord_length < rc:
            x1 += h
            y1 = self.at(x1)
            chord_length = ((x-x1)**2 + (y-y1)**2)**0.5

        if chord_length < rc:
            raise Exception("GertlerEnvelope: Unable to find the trailing edge intercept for the given parameters.")

        return x1, y1
    
    # Set the length of the envelope
    def set_length (self, l):
        self.length = l
        self.diameter = l * self.diameter / self.length

    # Returns the coefficients in an ordered way numpy usually accepts.
    @property
    def ordered_coeffs (self):
        return self.coeffs[::-1] + [0]

    # Returns the coefficients of the Gertler polynomial from standard parameters.
    def get_coefficients (params):
        (m, r0, r1, cp, _) = params

        A = np.array([   
            [1, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1],
            [m, m**2, m**3, m**4, m**5, m**6],
            [1, 2*m, 3*m**2, 4*m**3, 5*m**4, 6*m**5],
            [1, 2, 3, 4, 5, 6],
            [1/2, 1/3, 1/4, 1/5, 1/6, 1/7]
        ], float)

        B = np.array([2*r0, 0, 1/4, 0, -2*r1, 1/4*cp], float).T
        X = np.linalg.solve(A, B)
        
        return np.round(X, 4)
    
    # Returns a GertlerEnvelope from standard parameters.
    #
    # params = (m1, r0, r1, cp, l2d)
    # length = Length of the envelope
    # n = Number of points to be generated along the length
    def from_parameters (params, length, n):
        coeffs = GertlerEnvelope.get_coefficients(params)
        diameter = length / params[4]
        return GertlerEnvelope(coeffs, length, diameter, n)
    
    # Returns a GertlerEnvelope from standard parameters but with volume
    #
    # params = (m1, r0, r1, cp, l2d)
    # volume = Volume of the envelope
    # n = Number of points to be generated along the length
    # lobe_number = Number of lobes 
    # e, f, g = Multi lobe distances
    #
    # NOTE: The following code assumes all the lobes are of same shape and length.
    def from_parameters_volume (params, volume, n, lobe_number = 1, e = 0, f = 0, g = 0):
        coeffs = GertlerEnvelope.get_coefficients(params)
        cp = params[3]
        l2d = params[4]

        # Calculation of length using V = pi * D^2 * L * cp/4 for a monolobe
        length = ((4 * volume * l2d**2) / (np.pi * cp))**(1/3)

        if lobe_number == 1:
            # In case of monolobe, just pass the envelope
            return GertlerEnvelope(coeffs, length, length / l2d, n)
        else:
            # In case of multilobe, the formula for volume is complicated and scipy's fsolve is used to compute 
            # the required length for the given volume. 
            # 
            # For initial estimate of length, it is assumed there are no intersections between the lobe. This gives
            # us a good initial estimate making it more faster for us to arrive at the actual value.
            # 
            # Volume = Number of lobes * pi * D^2 * L_estimate * cp/4
            initial_estimate = length / (lobe_number**(1/3))
            initial_envelope = GertlerEnvelope(coeffs, initial_estimate, initial_estimate / l2d, n)

            # Different function is used in different case here. If we were to put this if condition within one function,
            # it will be computationally inefficient.
            if lobe_number == 2:
                # In case of bilobe design,
                def func (l):
                    # Update the new length to the envelope.
                    initial_envelope.set_length(l[0])
                    # New estimate of the volume from the new length.
                    volume_iter = initial_envelope.volume_bilobe(f)
                    # Returns a factor indicating how much did the volume change with the new length.
                    # With this factor approaching 0, the length will approach the required length accordingly.
                    return (volume - volume_iter) / volume  
            else:
                # In case of trilobe design,
                def func (l):
                    initial_envelope.set_length(l[0])
                    volume_iter = initial_envelope.volume_trilobe(e, f, g)
                    return (volume - volume_iter) / volume
                
            fsolve(func, initial_estimate)
            return initial_envelope

# Returns the half thickness of a symmetric NACA 4 digit airfoil.
#
# t = Maximum thickness as percentage of chord
def naca_airfoil_half_thickness_at (t, x):
    return 5 * t * (0.2969*x**0.5 - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)

# Returns an array of points representing a symmetric NACA airfoil.
#
# NOTE: These points are linearly scalable with length.
def naca_airfoil_points (t, n, l = 1):
    t *= l/100

    for i in range(n):
        x = i / n
        y = naca_airfoil_half_thickness_at(t, x)
        yield l * x, y
        
    yield l, 0

    # TODO: This is computing the same points for the lower surface again which maybe
    # inefficient computation wise. If we were to create an array and return, it would 
    # waste memory. Try to find a better solution for this.
    for i in range(1, n-1):
        x = 1 - i / n
        y = naca_airfoil_half_thickness_at(t, x)
        yield l * x, -y

    yield 0, 0

def length_scaling (volume_function, target_volume, initial_estimate):
    def f (l):
        volume = volume_function(l)
        return (target_volume - volume) / target_volume
    
    return fsolve(f, initial_estimate)[0]