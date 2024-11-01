from dataclasses import dataclass
from enum import Enum
from math import sqrt
from typing import Optional
import numpy as np
from sys import exit
from .standardsizes import Standard
from abc import ABC, abstractmethod

"""
MATERIALS
Unconfined Concrete     Create a material from Mander-Unconfined model given standard material properties.
Confined Concrete       Create a material from Mander-Confined model, given confinement pressure.
                            can be calculated with osfyber.conf_pressure()
Steel                   Create a material from E, Esh, ey, and other common steel properties
User Defined            Create an arbitrary material given several points of (strain,stress)
"""


class Color(Enum):
    Gray = 0
    Purple = 1
    Pink = 2
    Blue = 3
    Black = 4
    Green = 5
    Yellow = 6
    Red = 7
    White = 8
    Orange = 9


class Material(ABC):
    def __init__(self) -> None:
        self.fail = ""
        self.state: Color = Color.White
        self.useful_points: list[float] = []

    @abstractmethod
    def stress(self, ec: float) -> float:
        pass

    def stress_state(self, ec: float) -> tuple[float, Color]:
        stress = self.stress(ec)
        return (stress, self.state)

    @abstractmethod
    def tension(self) -> float:
        pass


@dataclass
class ConfinedPressureCircleProps:
    fyh: float
    bar_number: int
    D: float
    s: float


def conf_pressure_circle(props: ConfinedPressureCircleProps) -> float:
    """
    Confining pressure from a circular section
    """
    bar = Standard().bar(props.bar_number)
    ke = (1 - props.s / (2 * props.D)) ** 2
    fple = ke * 2 * bar.area * props.fyh / (props.D * props.s)
    return fple


@dataclass
class ConfinedPressureRectProps:
    fyh: float
    bar_number: int
    s: float
    b: float
    w: float
    nx: int
    ny: int


def conf_pressure_rect(props: ConfinedPressureRectProps) -> float:
    """
    Confining pressure from a rectangular section
    """
    # Is s_prime suppose to be hard coded?
    s_prime = 2
    bar = Standard().bar(props.bar_number)
    hx = props.b
    hy = props.w
    a_shx = bar.area * props.nx
    a_shy = bar.area * props.ny
    sum_w_sqr = (
        2 * (props.nx - 1) * (hx / (props.nx - 1) - bar.db) ** 2
        + 2 * (props.ny - 1) * (hy / (props.ny - 1) - bar.db) ** 2
    )
    # CONFINEMENT PROPERTIES
    ke = (
        (1 - (sum_w_sqr / (6 * hx * hy)))
        * (1 - s_prime / (2 * hx))
        * (1 - props.s / (2 * hy))
    )
    fpl_x = a_shx * props.fyh / (hy * props.s)
    fpl_y = a_shy * props.fyh / (hx * props.s)
    fple = ke * max(sqrt(fpl_x * fpl_y), 0.25 * max(fpl_x, fpl_x))
    # TODO Not sure what these are for
    # ke = (1 - s/(2*db) )**2
    # fple = ke*2*bar.area*fyh/(db*s)
    return fple


class ReinforcementProperties:
    """Store a list of bar locations as points"""

    def __init__(self, points: list[tuple[float, float]], bar_num: int, mat_id: int):
        bar = Standard().bar(bar_num)
        self.radius = bar.db / 2
        self.area = bar.area
        self.points = points
        self.mat_id = mat_id


@dataclass
class UnconfinedConcreteProps:
    fpc: float
    Ec: Optional[float] = None
    ecp: Optional[float] = None
    ecu: Optional[float] = None
    use_tension: bool = False


class UnconfConcMat(Material):
    """Model behavior from given concrete material properties"""

    def __init__(self, props: UnconfinedConcreteProps):
        self.fail = ""
        self.state: Color = Color.White
        self.fpc = props.fpc
        # General concrete properties
        # 1 psi == 0.006894757 MPa
        ratio = 1000
        # Testing metric input conversion
        # ratio = 0.006894757
        if props.Ec:
            self.Ec = props.Ec
        else:
            self.Ec = sqrt(self.fpc * ratio) * 57
        self.fcr = 4 * sqrt(self.fpc * ratio) / ratio

        # Unconfined properties
        # 3600 psi == 24.8211 MPa
        m = 1 + 3600 / (self.fpc * ratio) / ratio
        self.ecp = props.ecp if props.ecp else -self.fpc * m / self.Ec
        # self.r_unconf = 1/(1+self.fpc/(self.Ec*self.ecp))
        self.r_unconf = 2

        self.ecu = props.ecu if props.ecu else -0.005

        self.use_tension = props.use_tension

        # Important points to plot (In between colors) - also gives min/max strain
        self.useful_points = [self.ecu, self.ecp * 2, self.ecp, 0]

    def stress(self, ec: float) -> float:
        """Output unconfined concrete stress for a given strain"""
        self.fail = ""
        if ec > 0:
            # Tension Branch
            if self.use_tension:
                return ec  # not sure if this is right
            else:
                self.state = Color.White
                return 0
        else:
            # Compression Branch
            if ec >= self.ecp * 2:
                # Below ultimate unconfined strain
                r = self.r_unconf
                ecp = self.ecp
                fc = -r * (ec / ecp) / (r - 1 + (ec / ecp) ** r) * self.fpc
                if ec > self.ecp:
                    self.state = Color.Blue
                else:
                    self.state = Color.Purple
                return fc
            elif ec >= self.ecu:
                slope = self.stress(2 * self.ecp) / (2 * self.ecp - self.ecu)
                b = -slope * self.ecu
                self.state = Color.Red
                return slope * ec + b
            else:
                self.state = Color.White
                # Don't need to fail if unconf fails
                return 0

    def tension(self) -> float:
        """Need to check if user wants tension, if they do then use something like:
        elastic = ec * self.Ec
        if elastic > self.fcr:
            if ec <= 0.002:
                return 0.7*self.fcr/(1+sqrt(500*ec))
            else:
                # Ultimate tensile
                return 0
        else:
            return elastic
        """
        self.state = Color.White
        return 0


@dataclass
class ConfinedConcreteProps:
    fpc: float
    fple: float
    Ec: Optional[float] = None
    epcc: Optional[float] = None
    use_tension: bool = False


class ConfConcMat(Material):
    """Model behavior from given concrete material properties"""

    def __init__(self, props: ConfinedConcreteProps):
        self.fail = ""
        self.state: Color = Color.White
        self.fpc = props.fpc
        self.fple = props.fple
        # General concrete properties
        # 1 psi == 0.006894757 MPa
        ratio = 1000
        if props.Ec:
            self.Ec = props.Ec
        else:
            self.Ec = sqrt(self.fpc * ratio) * 57
        self.fcr = 4 * sqrt(self.fpc * ratio) / ratio

        # Confined properties
        # k1 = 5.4/sqrt(1+5*(self.fple/self.fpc))
        # self.fpcc = -self.fpc*(1+k1*(self.fple/self.fpc))
        # self.fpcc = -self.fpc*ratio * (2.254*sqrt(1 + (7.94*self.fple)/self.fpc)-(2*self.fple / self.fpc) - 1.254)
        self.fpcc = -6.899

        # TODO HARDCODED ASSUME 0.002????
        # self.epc0 = -0.002219
        self.epc0 = -0.002

        # m = 1+3600/(self.fpc*ratio)/ratio
        if props.epcc:
            self.epcc = props.epcc
        else:
            self.epcc = self.epc0 * (1 + 5 * (-self.fpcc / self.fpc - 1))

        self.Esec = self.fpcc / self.epcc
        self.r_conf = self.Ec / (self.Ec - self.Esec)
        # self.ecu = self.epcc*(1+20*(self.fple/self.fpc))
        self.ecu = -(0.004 + 1.4 * 0.00831 * 68 * 0.09 / 6.899)
        # Debug Variables
        # print(self.__dict__)

        self.use_tension = props.use_tension

        # Important points to plot (In between colors) - also gives min/max strain
        self.useful_points = [self.ecu, self.epcc, 0, self.fcr, 0.002]

    def stress(self, ec: float) -> float:
        """Output confined concrete stress for a given strain"""
        self.fail = ""
        if ec > 0:
            self.state = Color.White
            return 0
            # TODO Tension Branch - might need this if user specifies
            """
            self.state = Color.Blue
            elastic = ec * self.Ec
            if elastic > self.fcr:
                if ec <= 0.002:
                    return 0.7*self.fcr/(1+sqrt(500*ec))
                else:
                    # Ultimate tensile
                    self.state = Color.White
                    return 0
            else:
                return elastic
            """
        else:
            # Compression Branch
            if ec >= self.ecu:
                r = self.r_conf
                x = ec / self.epcc
                # FIXME - This is a temporary fix for the division by zero error
                fcc: float = self.fpcc * x * r / (r - 1 + x**r)
                if ec > self.epcc:
                    self.state = Color.Green
                else:
                    self.state = Color.Orange
                return fcc
            else:
                # Over ultimate strain - Hoop fracture
                self.state = Color.Black
                self.fail = (
                    f"Confined Concrete Crushing\n"
                    f"\tMax Available Strain={round(self.ecu,5)}\n"
                    f"\tStrain Experienced={round(ec,5)}"
                )
                return 0

    def tension(self) -> float:
        self.state = Color.White
        return 0


@dataclass
class SteelProps:
    E: float
    fy: float
    fsu: float
    e_sh: float
    e_su: float
    P: float


class SteelMat(Material):
    """Model behavior from given steel material properties"""

    def __init__(self, props: SteelProps):
        self.fail = ""
        self.state: Color = Color.Black
        # Needed steel properties
        self.Es = props.E
        self.fy = props.fy
        self.fsu = props.fsu
        self.e_y = props.fy / props.E
        self.e_sh = props.e_sh
        self.e_su = props.e_su
        self.P = props.P

        # Important points to plot (In between colors) - also gives min/max strain
        self.useful_points = [
            -self.e_su,
            -self.e_sh,
            -self.e_y,
            0,
            self.e_y,
            self.e_sh,
            self.e_su,
        ]

    def stress(self, e: float) -> float:
        """Output stress in steel from a given strain"""
        # TODO - Combine ten/comp, multiply by strain/abs(strain) for sign
        self.fail = ""
        if e < 0:
            # Mirrored Response of tension
            if e > -self.e_y:
                # Elastic
                self.state = Color.Yellow
                return self.Es * e
            elif e > -self.e_sh:
                # Pseudo Plastic Before SH
                self.state = Color.Pink
                return -self.fy
            elif e > -self.e_su:
                # Strain Hardening
                fsu = self.fsu
                e_su = self.e_su
                self.state = Color.Red
                # FIXME - This is a temporary fix for complex numbers
                sh1: float = ((e_su + e) / (e_su - self.e_sh)) ** self.P
                return -(fsu - (fsu - self.fy) * sh1)
            else:
                # Tension Fracture
                self.state = Color.Black
                self.fail = (
                    f"Steel Fracture\n"
                    f"\tMax Available Strain={abs(round(self.e_su,5))}\n"
                    f"\tStrain Experienced={abs(round(e,5))}"
                )
                return 0
        else:
            # Tension Branch
            if e < self.e_y:
                # Elastic
                self.state = Color.Yellow
                return self.Es * e
            elif e < self.e_sh:
                # Pseudo Plastic Before SH
                self.state = Color.Pink
                return self.fy
            elif e < self.e_su:
                # Strain Hardening
                fsu = self.fsu
                e_su = self.e_su
                self.state = Color.Red
                # FIXME - This is a temporary fix for complex numbers
                sh2: float = ((e_su - e) / (e_su - self.e_sh)) ** self.P
                return fsu - (fsu - self.fy) * sh2
            else:
                # Tension Fracture
                self.state = Color.Black
                self.fail = (
                    f"Steel Fracture\n"
                    f"\tMax Available Strain={abs(round(self.e_su,5))}\n"
                    f"\tStrain Experienced={abs(round(e,5))}"
                )
                return 0

    def tension(self) -> float:
        self.state = Color.White
        return 0


@dataclass
class UserMaterialProps:
    points: list[tuple[float, float]]
    mirror: bool = False


class UserMat(Material):
    """Model behavior from given (strain,stress) points"""

    def __init__(self, props: UserMaterialProps):
        self.fail = ""
        self.state: Color = Color.White

        self.strains = [i[0] for i in props.points]
        self.stresses = [i[1] for i in props.points]
        self.yield_strain = abs(self.strains[1])
        self.mirror = props.mirror
        # flip if given in reverse order (needed for material plots)
        if self.strains[0] == max(self.strains):
            self.strains = self.strains[::-1]
            self.stresses = self.stresses[::-1]
        if props.mirror:
            # TODO: user input for failure levels (min/max/both)
            # copy points over line y=-x
            if min(self.strains) < 0:
                print(
                    "You can't mirror a user material if there are negative strain values"
                )
                exit()
            # if there's a point of 0 strain, don't copy it
            # or maybe there's a np function to remove duplicates?
            if min(self.strains) == 0:
                self.strains = [-i for i in self.strains[::-1]] + self.strains[1:]
                self.stresses = [-i for i in self.stresses[::-1]] + self.stresses[1:]
            else:
                self.strains = [-i for i in self.strains[::-1]] + self.strains
                self.stresses = [-i for i in self.stresses[::-1]] + self.stresses

        # Important points to plot (In between colors) - also gives min/max strain
        self.useful_points = self.strains

    def stress(self, e: float) -> float:
        """Output stress from a given strain"""
        self.fail = ""
        if (e < min(self.strains)) or (e > max(self.strains)):
            # No more strength - below min or above max strain
            if self.mirror:
                # stop analysis if user material is mirrored
                self.fail = (
                    f"User Material Strain Limit\n"
                    f"\tMax Available Strain={abs(round(self.strains[0],5))}\n"
                    f"\tStrain Experienced={abs(round(e,5))}"
                )
            self.state = Color.Black
            return 0
        if abs(e) < self.yield_strain:
            self.state = Color.Yellow
        else:
            self.state = Color.Red
        # Todo: Want to color with gradient?
        # if self.gradient:
        #    loc = abs(e) / max(max(self.strains),-min(self.strains))
        #    self.state = (loc,0,1-loc)
        v = np.interp(e, self.strains, self.stresses)
        return float(v[0])

    def tension(self) -> float:
        self.state = Color.White
        return 0
