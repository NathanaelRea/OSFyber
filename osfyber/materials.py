from math import sqrt
import numpy as np
from sys import exit
from osfyber.standardsizes import Standard

"""
MATERIALS
Unconfined Concrete     Create a material from Mander-Unconfined model given standard material properties.
Confined Concrete       Create a material from Mander-Confined model, given confinement pressure.
                            can be calculated with osfyber.conf_pressure()
Steel                   Create a material from E, Esh, ey, and other common steel properties
User Defined            Create an arbitrary material given several points of (strain,stress)
"""


def conf_pressure_circle(fyh, bar_number, dia, s):
    """
    Confining pressure from a circular section
    """
    bar = Standard().bar(bar_number)
    ke = (1 - s/(2*dia))**2
    fple = ke * 2 * bar.area * fyh/(dia*s)
    return fple


def conf_pressure_rect(fyh, bar_number, s, b, w, nx, ny):
    """
    Confining pressure from a rectangular section
    """
    # Is s_prime suppose to be hard coded?
    s_prime = 2
    bar = Standard().bar(bar_number)
    hx = b
    hy = w
    a_shx = bar.area * nx
    a_shy = bar.area * ny
    sum_w_sqr = 2*(nx-1)*(hx/(nx-1)-bar.db)**2 + 2*(ny-1)*(hy/(ny-1)-bar.db)**2
    # CONFINEMENT PROPERTIES
    ke = (1-(sum_w_sqr/(6*hx*hy)))*(1-s_prime/(2*hx))*(1-s/(2*hy))
    fpl_x = a_shx*fyh / (hy*s)
    fpl_y = a_shy*fyh / (hx*s)
    fple = ke * max(sqrt(fpl_x*fpl_y), 0.25*max(fpl_x, fpl_x))
    # TODO Not sure what these are for
    # ke = (1 - s/(2*db) )**2
    # fple = ke*2*bar.area*fyh/(db*s)
    return fple


class ReinforcementProperties:
    """Store a list of bar locations as points"""
    def __init__(self, points, bar_num, mat_id):
        bar = Standard().bar(bar_num)
        self.radius = bar.db/2
        self.area = bar.area
        self.points = points
        self.mat_id = mat_id


class UnconfConcMat:
    """Model behavior from given concrete material properties"""
    def __init__(self, kwargs):
        self.fail = False
        self.state = 'White'
        self.fpc = kwargs['fpc']
        # General concrete properties
        # 1 psi == 0.006894757 MPa
        ratio = 1000
        # Testing metric input conversion
        # ratio = 0.006894757
        if 'Ec' in kwargs:
            self.Ec = kwargs['Ec']
        else:
            self.Ec = sqrt(self.fpc * ratio) * 57
        self.fcr = 4*sqrt(self.fpc*ratio)/ratio
        
        # Unconfined properties
        # 3600 psi == 24.8211 MPa
        m = 1+3600/(self.fpc*ratio)/ratio
        if 'ecp' in kwargs:
            self.ecp = kwargs['ecp']
        else:
            self.ecp = -self.fpc*m/self.Ec
        # self.r_unconf = 1/(1+self.fpc/(self.Ec*self.ecp))
        self.r_unconf = 2
        
        if 'ecu' in kwargs:
            self.ecu = kwargs['ecu']
        else:
            self.ecu = -0.005
        
        if 'tension' in kwargs:
            self.tension = kwargs['tension']
        else:
            self.tension = False
        
        # Important points to plot (In between colors) - also gives min/max strain
        self.useful_points = [self.ecu, self.ecp*2, self.ecp, 0]
        
    def stress(self, ec):
        """Output unconfined concrete stress for a given strain"""
        self.fail = False
        if ec > 0:
            # Tension Branch
            if self.tension:
                return ec  # not sure if this is right
            else:
                self.state = 'White'
                return 0
        else:
            # Compression Branch
            if ec >= self.ecp * 2:
                # Below ultimate unconfined strain
                r = self.r_unconf
                ecp = self.ecp
                fc = -r*(ec/ecp)/(r-1 + (ec/ecp)**r)*self.fpc
                if ec > self.ecp:
                    self.state = 3
                else:
                    self.state = 1
                return fc
            elif ec >= self.ecu:
                slope = self.stress(2*self.ecp)/(2*self.ecp - self.ecu)
                b = -slope*self.ecu
                self.state = 7
                return slope*ec+b
            else:
                self.state = 8
                # Don't need to fail if unconf fails
                return 0
            
    # def tension(self, ec):
    def tension(self):
        """Output tension concrete stress for a given strain"""
        """ Need to check if user wants tension, if they do then use something like:
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
        self.state = 8
        return 0


class ConfConcMat:
    """Model behavior from given concrete material properties"""
    def __init__(self, kwargs):
        self.fail = False
        self.state = 'White'
        self.fpc = kwargs['fpc']
        self.fple = kwargs['fple']
        # General concrete properties
        # 1 psi == 0.006894757 MPa
        ratio = 1000
        if 'Ec' in kwargs:
            self.Ec = kwargs['Ec']
        else:
            self.Ec = sqrt(self.fpc * ratio) * 57
        self.fcr = 4*sqrt(self.fpc*ratio)/ratio
        
        # Confined properties
        # k1 = 5.4/sqrt(1+5*(self.fple/self.fpc))
        # self.fpcc = -self.fpc*(1+k1*(self.fple/self.fpc))
        # self.fpcc = -self.fpc*ratio * (2.254*sqrt(1 + (7.94*self.fple)/self.fpc)-(2*self.fple / self.fpc) - 1.254)
        self.fpcc = -6.899
        
        # TODO HARDCODED ASSUME 0.002????
        # self.epc0 = -0.002219
        self.epc0 = -0.002
        
        # m = 1+3600/(self.fpc*ratio)/ratio
        if 'epcc' in kwargs:
            self.epcc = kwargs['ecp']
        else:
            self.epcc = self.epc0*(1+5*(-self.fpcc/self.fpc - 1))
        
        self.Esec = self.fpcc / self.epcc
        self.r_conf = self.Ec/(self.Ec - self.Esec)
        # self.ecu = self.epcc*(1+20*(self.fple/self.fpc))
        self.ecu = - (0.004 + 1.4 * 0.00831 * 68 * 0.09 / 6.899)
        # Debug Variables
        # print(self.__dict__)
        
        if 'tension' in kwargs:
            self.tension = kwargs['tension']
        else:
            self.tension = False
        
        # Important points to plot (In between colors) - also gives min/max strain
        self.useful_points = [self.ecu, self.epcc, 0, self.fcr, 0.002]
        
    def stress(self, ec):
        """Output confined concrete stress for a given strain"""
        self.fail = False
        if ec > 0:
            self.state = "White"
            return 0
            # TODO Tension Branch - might need this if user specifies
            # noinspection PyUnreachableCode
            """
            self.state = 'Cyan'
            elastic = ec * self.Ec
            if elastic > self.fcr:
                if ec <= 0.002:
                    return 0.7*self.fcr/(1+sqrt(500*ec))
                else:
                    # Ultimate tensile
                    self.state = 8
                    return 0
            else:
                return elastic
            """
        else:
            # Compression Branch
            if ec >= self.ecu:
                r = self.r_conf
                x = ec / self.epcc
                fcc = self.fpcc*x*r/(r-1+x**r)
                if ec > self.epcc:
                    self.state = 'Green'
                else:
                    self.state = 'Orange'
                return fcc
            else:
                # Over ultimate strain - Hoop fracture
                self.state = 'Black'
                self.fail = f"Confined Concrete Crushing\n" \
                            f"\tMax Available Strain={round(self.ecu,5)}\n" \
                            f"\tStrain Experienced={round(ec,5)}"
                return 0


# noinspection PyPep8Naming
class SteelMat:
    """Model behavior from given steel material properties"""
    def __init__(self, E, fy, fsu, e_sh, e_su, P):
        self.fail = False
        self.state = 'Black'
        # Needed steel properties
        self.Es = E
        self.fy = fy
        self.fsu = fsu
        self.e_y = fy/E
        self.e_sh = e_sh
        self.e_su = e_su
        self.P = P
        
        # Important points to plot (In between colors) - also gives min/max strain
        self.useful_points = [-self.e_su, -self.e_sh, -self.e_y, 0, self.e_y, self.e_sh, self.e_su]
        
    def stress(self, e):
        """Output stress in steel from a given strain"""
        # TODO - Combine ten/comp, multiply by strain/abs(strain) for sign
        self.fail = False
        if e < 0:
            # Mirrored Response of tension
            if e > -self.e_y:
                # Elastic
                self.state = 'Yellow'
                return self.Es * e
            elif e > -self.e_sh:
                # Pseudo Plastic Before SH
                self.state = 'Pink'
                return -self.fy
            elif e > -self.e_su:
                # Strain Hardening
                fsu = self.fsu
                e_su = self.e_su
                self.state = 'Red'
                return -(fsu-(fsu-self.fy)*((e_su+e)/(e_su-self.e_sh))**self.P)
            else:
                # Tension Fracture
                self.state = 'Black'
                self.fail = f"Steel Fracture\n" \
                            f"\tMax Available Strain={abs(round(self.e_su,5))}\n" \
                            f"\tStrain Experienced={abs(round(e,5))}"
                return 0
        else:
            # Tension Branch
            if e < self.e_y:
                # Elastic
                self.state = 'Yellow'
                return self.Es * e
            elif e < self.e_sh:
                # Pseudo Plastic Before SH
                self.state = 'Pink'
                return self.fy
            elif e < self.e_su:
                # Strain Hardening
                fsu = self.fsu
                e_su = self.e_su
                self.state = 'Red'
                return fsu-(fsu-self.fy)*((e_su-e)/(e_su-self.e_sh))**self.P
            else:
                # Tension Fracture
                self.state = 'Black'
                self.fail = f"Steel Fracture\n" \
                            f"\tMax Available Strain={abs(round(self.e_su,5))}\n" \
                            f"\tStrain Experienced={abs(round(e,5))}"
                return 0


class UserMat:
    """Model behavior from given (strain,stress) points"""
    def __init__(self, points, mirror=False):
        self.fail = False
        self.state = 'White'

        self.strains = [i[0] for i in points]
        self.stresses = [i[1] for i in points]
        self.yield_strain = abs(self.strains[1])
        self.mirror = mirror
        # flip if given in reverse order (needed for material plots)
        if self.strains[0] == max(self.strains):
            self.strains = self.strains[::-1]
            self.stresses = self.stresses[::-1]
        if mirror:
            # TODO: user input for failure levels (min/max/both)
            # copy points over line y=-x
            if min(self.strains) < 0:
                print("You can't mirror a user material if there are negative strain values")
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
            
    def stress(self, e):
        """Output stress from a given strain"""
        self.fail = False
        if (e < min(self.strains)) or (e > max(self.strains)):
            # No more strength - below min or above max strain
            if self.mirror:
                # stop analysis if user material is mirrored
                self.fail = f"User Material Strain Limit\n" \
                            f"\tMax Available Strain={abs(round(self.strains[0],5))}\n" \
                            f"\tStrain Experienced={abs(round(e,5))}"
            self.state = 'Black'
            return 0
        if abs(e) < self.yield_strain:
            self.state = 6
        else:
            self.state = 7
        # Todo: Want to color with gradient?
        # if self.gradient:
        #    loc = abs(e) / max(max(self.strains),-min(self.strains))
        #    self.state = (loc,0,1-loc)
        return np.interp(e, self.strains, self.stresses)
