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

import numpy as np

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
    
    # Iterator to get the radial points along the length of envelope.
    def points (self, truncation = 0):
        X = np.linspace(0, 1 - truncation, self.n)
        R = np.polyval(self.coeffs[::-1] + [0], X)
        R[R < 0] = 0

        # In case if there is no truncation, the final point must lie on the axis.
        if not truncation:
            R[-1] = 0 

        return zip(X * self.length, self.diameter * np.sqrt(R))

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

        return (x1, y1)
    
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