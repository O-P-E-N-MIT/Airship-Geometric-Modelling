# Input Variables for Gertler Envelope Generation
#
# m1 = Point of maximum thickness
# r0 = Nose radius at the front of the airship
# r1 = Tail radius at the back of the airship
# cp = Prismatic coefficient
# l2d = Length to diameter ratio

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

    def __init__ (self, params, n, l = 1):
        self.params = params
        self.coeffs = GertlerEnvelope.get_coefficients(params)
        self.n = n
        self.length = l
        self.diameter = l / params[4]

    def at (self, x):
        x = x / self.length
        y = self.diameter * (self.coeffs[0]*x + self.coeffs[1]*x**2 + self.coeffs[2]*x**3 + self.coeffs[3]*x**4 + self.coeffs[4]*x**5 + self.coeffs[5]*x**6)**0.5
        return y

    def points (self):
        for i in range(self.n):
            x = i / self.n
            y = self.diameter * (self.coeffs[0]*x + self.coeffs[1]*x**2 + self.coeffs[2]*x**3 + self.coeffs[3]*x**4 + self.coeffs[4]*x**5 + self.coeffs[5]*x**6)**0.5
            yield (self.length * x, y)

        yield (self.length, 0)

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

        return (x1, y1)

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

        return X

# Input Variables for Symmetric NACA Airfoil
#
# t = Maximum thickness as percentage of chord

def naca_airfoil_camber (t, x):
    return 5 * t * (0.2969*x**0.5 - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)

# NOTE: These points are linearly scalable with length.
# 
# I did not make this an iterator function like the Gertler Envelope because 
# the same points can be used to scale for the tapering wing.
def naca_airfoil_points (t, n, l = 1):
    t *= l/100
    points = []

    for i in range(n):
        x = i / n
        y = naca_airfoil_camber(t, x)
        points.append((l * x, y))
    
    points.append((l, 0))

    n = len(points)
    for i in range(1, n-1):
        points.append((points[n-(i + 1)][0], -points[n-(i + 1)][1]))

    points.append(points[0])

    return points