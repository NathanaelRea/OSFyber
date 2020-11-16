import sys


def secant_method(fct, _value, tol=10 ** -2, _iter=10 ** -10):
    """
    Minimize a given function to zero with the Secant method
    I'm not sure how to find the derivative of a black box function,
    And so we can't directly find the derivative, we must first find
    The linear slope between two evaluation points
    fct:    an actual python function with one input/output 
        This is set up in advance to minimize to zero
    start:    the first input
    tol:    The tolerance by witch to quit when the function return
        is < |fct ret|
    iter:    How much to step by in each loop
    """
    # TODO remove secant and bisection method functions, use numpy
    initial_val = _value
    for step in range(10):
        # Get current location within function
        cur = fct(_value)
        # We see if we have minimized the function to zero
        if abs(cur) < tol:
            return _value
        _value += _iter
        slope = (fct(_value) - cur) / _iter
        # If slope == 0 -> Divergence
        # We did not pick an initial value within the zone of influence
        if slope == 0:
            return bisection_method(fct, initial_val, tol)
        _value -= cur / slope
    # Secant method should be able to converge within 4 or so steps.
    return bisection_method(fct, initial_val, tol)


def bisection_method(fct, i_val, tol):
    """
    Sometimes, the secant/Newton Rhapshon method can diverge
    This happens when the initial guess is not within the zone of influence
    Bisection is guaranteed to always converge, since we start with two
    sides of the function that return pos and neg values
    """
    # Find a scope - Function values on either side have opposite signs
    scope = 5
    while True:
        # Check if we get +/- on either side of the scope
        a = i_val - scope
        b = i_val + scope
        left = fct(a)
        right = fct(b)
        # check if opposite signs
        if abs(left) / left != abs(right) / right:
            break
        scope *= 2
        if scope > 10 ** 3:
            print("Bisection Method Scope out of bounds\nExternally applied force may be larger than section capacity?")
            sys.exit()
    # Try 100 iterations
    a_val = fct(a)
    c_val = fct((a + b) / 2)
    n_max = 100
    for n in range(1, n_max):  # limit iterations to prevent infinite loop
        c = (a + b) / 2  # new midpoint
        if abs(c_val - a_val) / 2 <= tol:  # solution found
            return c
        a_val = fct(a)
        c_val = fct(c)
        # Check if midpoint and left are same sign
        if c_val / abs(c_val) == a_val / abs(a_val):
            a = c
        else:
            b = c  # new interval
    print("Bisection Method Failed")
    sys.exit()


class State:
    """State information for a phi step"""

    def __init__(self, mat_state, max_strain, min_strain, y_loc):
        self.mat_state = mat_state  # List of state values for each patch
        self.max_strain = max_strain
        self.min_strain = min_strain
        self.strains = []  # List of strain values for each patch
        self.stresses = []  # List of strain values for each patch
        self.y_loc = y_loc
        # TODO REMOVE MAX_STRAIN/MIN_STRAIN replace with max(self.strains)


class Fiber:
    """
    A fiber is just a small piece of geometry.
    It's only properties are area, cords, and material.
    """

    def __init__(self, area, cords, mat_id):
        self.area = area
        self.xy = cords
        self.mat_id = mat_id
        self.fail = False
        # self.distance = None? # for quick calculations


class FiberModel:
    """
    A method to put all mesh and reinforcement information in one place
    That way we can easily calculate the state of the system, and results
    Such as M, and force balance at a given curvature level
    """

    def __init__(self):
        """Analysis Model Class
        Functions to find equilibrium and return state data/ M from phi level"""
        self.fibers = {}
        self.materials = {}
        # location where strain goes from negative to positive
        self.zero_strain_location = 0
        self.phi = 0
        # TODO CAN I JUST FIND THE DISTANCE HERE (once) from (0,0) and given angle?
        self.P = 0
        self.fail = False
        self.maxy = 0
        self.miny = 0
        # State Object, used for max/min strains of this particular fiber
        self.state = None

    def calc_strain(self, cords, y_intercept):
        # TODO - ADD CURVATURE ANGLE from +y
        # (Currently only does +y direction)
        # strain = self.phi * (y_intercept - cords[1])
        return self.phi * (y_intercept - cords[1])

    def force_balance(self, y_intercept):
        """
        Given the intercept (zero_strain_location) and slope (curvature)
        we can find the forces of each fiber and sum them (to be able to balance)
        """
        self.zero_strain_location = y_intercept
        sum_force = 0
        for i, fiber in self.fibers.items():
            strain = self.calc_strain(fiber.xy, y_intercept)
            stress = self.materials[fiber.mat_id].stress(strain)
            sum_force += fiber.area * stress
        # print(f"{y_intercept}: {sum_F}")
        return sum_force + self.P

    def calc_moment(self):
        """Find the internal moment of the centroid to balance the
        External moment and the internal stress from curvature"""
        # TODO: Change this to save state. (Or have another fct that calls this fct)
        # M is just one property of this solved curvature level

        # Go through all elements (mesh/reinforcement) and calculate strain from linear centroid, then stress
        # TODO - Find some transformation value c1*x+c2*y -> distance from strain reversal/centroid
        # BUT ONLY CALCULATE THAT ONCE>>>> THEN JUST CALCULATE THE Sin(theta) OR WHATEVER TO FIND THE STEP
        # WE DO NOT NEED TO RECALCULATE THIS EVERY TIME
        # BECAUSE WE HAVE A KNOWN STEP SIZE, AND MOMENT ROTATION DOES NOT CHANGE

        # Calculate strain at furthest fibers
        max_strain = self.phi * (self.zero_strain_location - self.maxy)
        min_strain = self.phi * (self.zero_strain_location - self.miny)

        sum_moment = 0
        mat_states = []
        strains = []
        stresses = []
        for i, fiber in self.fibers.items():
            strain = self.calc_strain(fiber.xy, self.zero_strain_location)
            stress = self.materials[fiber.mat_id].stress(strain)
            strains.append(strain)
            stresses.append(stress)
            sum_moment += fiber.area * stress * (self.zero_strain_location - fiber.xy[1])
            mat_states.append(self.materials[fiber.mat_id].state)
            fiber.fail = self.materials[fiber.mat_id].fail
            if fiber.fail:
                rnd_loc = [round(i, 3) for i in fiber.xy]
                self.fail = fiber.fail + f"\n\tMat_id={fiber.mat_id}\n\tLocation={rnd_loc}"
        # Add load P into moment calculation
        sum_moment += self.P * self.zero_strain_location
        self.state = State(mat_states, max_strain, min_strain, self.zero_strain_location)
        self.state.strains = strains
        self.state.stresses = stresses
        return sum_moment
