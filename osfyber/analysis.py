from typing import Optional
from osfyber.materials import Color, Material


class State:
    """State information for a phi step"""

    def __init__(
        self, mat_state: list[Color], max_strain: float, min_strain: float, y_loc: float
    ):
        self.mat_state: list[Color] = mat_state  # List of state values for each patch
        self.max_strain = max_strain
        self.min_strain = min_strain
        self.strains: list[float] = []  # List of strain values for each patch
        self.stresses: list[float] = []  # List of strain values for each patch
        self.y_loc = y_loc
        # TODO REMOVE MAX_STRAIN/MIN_STRAIN replace with max(self.strains)


class Fiber:
    """
    A fiber is just a small piece of geometry.
    It's only properties are area, cords, and material.
    """

    def __init__(self, area: float, cords: tuple[float, float], mat_id: int):
        self.area = area
        self.xy = cords
        self.mat_id = mat_id
        self.fail = ""
        # self.distance = None? # for quick calculations


class FiberModel:
    """
    A method to put all mesh and reinforcement information in one place
    That way we can easily calculate the state of the system, and results
    Such as M, and force balance at a given curvature level
    """

    def __init__(self) -> None:
        """Analysis Model Class
        Functions to find equilibrium and return state data/ M from phi level"""
        self.fibers: dict[int, Fiber] = {}
        self.materials: dict[int, Material] = {}
        # location where strain goes from negative to positive
        self.zero_strain_location = 0.0
        self.phi = 0.0
        # TODO CAN I JUST FIND THE DISTANCE HERE (once) from (0,0) and given angle?
        self.P = 0.0
        self.fail = ""
        self.max_y = 0.0
        self.min_y = 0.0
        # State Object, used for max/min strains of this particular fiber
        self.state: Optional[State] = None

    def calc_strain(self, cords: tuple[float, float], y_intercept: float) -> float:
        # TODO - ADD CURVATURE ANGLE from +y
        # (Currently only does +y direction)
        # strain = self.phi * (y_intercept - cords[1])
        return self.phi * (y_intercept - cords[1])

    def force_balance(self, y_intercept: float) -> float:
        """
        Given the intercept (zero_strain_location) and slope (curvature)
        we can find the forces of each fiber and sum them (to be able to balance)
        """
        self.zero_strain_location = y_intercept
        sum_force = 0.0
        for i, fiber in self.fibers.items():
            strain = self.calc_strain(fiber.xy, y_intercept)
            stress = self.materials[fiber.mat_id].stress(strain)
            sum_force += fiber.area * stress
        return sum_force + self.P

    def calc_moment(self) -> float:
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
        max_strain = self.phi * (self.zero_strain_location - self.max_y)
        min_strain = self.phi * (self.zero_strain_location - self.min_y)

        sum_moment = 0.0
        mat_states = []
        strains = []
        stresses = []
        for _, fiber in self.fibers.items():
            strain = self.calc_strain(fiber.xy, self.zero_strain_location)
            stress = self.materials[fiber.mat_id].stress(strain)
            strains.append(strain)
            stresses.append(stress)
            sum_moment += (
                fiber.area * stress * (self.zero_strain_location - fiber.xy[1])
            )
            mat_states.append(self.materials[fiber.mat_id].state)
            fiber.fail = self.materials[fiber.mat_id].fail
            if fiber.fail:
                rnd_loc = [round(i, 3) for i in fiber.xy]
                self.fail = (
                    fiber.fail + f"\n\tMat_id={fiber.mat_id}\n\tLocation={rnd_loc}"
                )
        # Add load P into moment calculation
        sum_moment += self.P * self.zero_strain_location
        self.state = State(
            mat_states, max_strain, min_strain, self.zero_strain_location
        )
        self.state.strains = strains
        self.state.stresses = stresses
        return sum_moment
