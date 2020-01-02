from math import sqrt
import numpy as np

"""
Materials and explanations:
Unconfined Concrete		Create a material from Mander-Unconfined model given standard material properties.
Confined Concrete		Create a material from Mander-Confined model, given confinement pressure.
							can be calculatled with osfyber.conf_pressure()
Steel					Create a material from E, Esh, ey, and other common steel properties
User Defined			Create an arbitrary material given several points of (strain,stress)
"""

def conf_pressure_circ(fyh, bar, dia, s):
	"""
	Confining pressure from a circular section
	"""
	bars = {3:(0.375,0.110), 4:(0.500,0.200), 5:(0.625,0.310), 
		6:(0.750,0.440), 7:(0.875,0.600), 8:(1.000,0.790), 9:(1.128,1.000), 
		10:(1.270,1.270), 11:(1.410,1.560), 14:(1.693,2.250), 18:(2.257,4.000)}
	Ab = bars[bar][1]
	ke = (1 - s/(2*dia) )**2
	fple = ke*2*Ab*fyh/(dia*s)
	return fple

def conf_pressure_rect(fyh, bar, s, b, w, nx, ny):
	"""
	Confining pressure from a rectangular section
	"""
	db = self.bars[bar][0]
	Ab = self.bars[bar][1]
	hx = b
	hy = w
	Ashx = Ab*nx
	Ashy = Ab*ny
	sum_w_sqr = 2*(nx-1)*(hx/(nx-1)-db)**2 + 2*(ny-1)*(hy/(ny-1)-db)**2
	### Confinement props
	ke   = (1-(sum_w_sqr/(6*hx*hy)))*(1-s_prime/(2*hx))*(1-s/(2*hy))
	fplx = Ashx*fyh / (hy*s)
	fply = Ashy*fyh / (hx*s)
	fple = ke * max(sqrt(fplx*fply), 0.25*max(fplx,fplx))
	return fplelf.bars[bar][1]
	ke = (1 - s/(2*dia) )**2
	fple = ke*2*Ab*fyh/(dia*s)
	return fple


class reinf_props:
	"""Store values for reinforcement ring/line"""
	def __init__(self, points, bar, mat_id):
		# Standard Bar Properties #:( db (in) ,Ab (in^2) )
		bars = {3:(0.375,0.110), 4:(0.500,0.200), 5:(0.625,0.310), 
			6:(0.750,0.440), 7:(0.875,0.600), 8:(1.000,0.790), 9:(1.128,1.000), 
			10:(1.270,1.270), 11:(1.410,1.560), 14:(1.693,2.250), 18:(2.257,4.000)}
		self.radius = bars[bar][0]/2
		self.area   = bars[bar][1]
		self.points = points
		self.mat_id = mat_id


class unconf_conc_material:
	"""Model behavior from given concrete material properties"""
	def __init__(self, kwargs):
		self.fpc = kwargs['fpc']
		### General concrete properties
		# 1 psi == 0.006894757 MPa
		ratio = 1000
		#ratio = 0.006894757
		if 'Ec' in kwargs:
			self.Ec = kwargs['Ec']
		else:
			self.Ec = sqrt(self.fpc * ratio) * 57
		self.fcr = 4*sqrt(self.fpc*ratio)/ratio
		
		### Unconf properties
		# 3600 psi == 24.8211 MPa
		m = 1+3600/(self.fpc*ratio)/ratio
		if 'ecp' in kwargs:
			self.ecp = kwargs['ecp']
		else:
			self.ecp = -self.fpc*m/self.Ec
		self.r_unconf = 1/(1+self.fpc/(self.Ec*self.ecp))
		
		if 'ecu' in kwargs:
			self.ecu = kwargs['ecu']
		else:
			self.ecu = self.ecp * 3
		
		self.min_strain = self.ecu
		self.max_strain = 0
		
	def stress(self, ec):
		"""Output unconfined concrete stress for a given strain"""
		self.fail = False
		if ec > 0:
			# Tension Branch
			return self.tension(ec)
		else:
			# Compression Branch
			if ec >= self.ecp *2:
				# Below ultimate unconfined strain
				r = self.r_unconf
				r = 2
				ecp = self.ecp
				fc = -r*(ec/ecp)/(r-1 + (ec/ecp)**r)*self.fpc
				if ec > self.ecp:
					self.state = 3
				else:
					self.state = 1
				return fc
			elif ec >=  self.ecu:
				slope = self.stress(2*self.ecp)/(2*self.ecp - self.ecu)
				b = -slope*self.ecu
				self.state = 7
				return slope*ec+b
			else:
				self.state = 8
				# Don't need to fail if unconf fails
				return 0
			
	def tension(self, ec):
		"""Output tension concrete stress for a given strain"""
		self.state = 8
		return 0
		"""
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

class conf_conc_material:
	"""Model behavior from given concrete material properties"""
	def __init__(self, kwargs):
		self.fpc = kwargs['fpc']
		self.fple = kwargs['fple']
		### General concrete properties
		# 1 psi == 0.006894757 MPa
		ratio = 1000
		if 'Ec' in kwargs:
			self.Ec = kwargs['Ec']
		else:
			self.Ec = sqrt(self.fpc * ratio) * 57
		self.fcr = 4*sqrt(self.fpc*ratio)/ratio
		
		### Confined properties
		#k1 = 5.4/sqrt(1+5*(self.fple/self.fpc))
		#self.fpcc = -self.fpc*(1+k1*(self.fple/self.fpc))
		#self.fpcc = -self.fpc*ratio * (2.254*sqrt(1 + (7.94*self.fple)/self.fpc)-(2*self.fple / self.fpc) - 1.254)
		self.fpcc = -6.899
		
		# TODO HARDCODED ASSUME 0.002????
		#self.epc0 = -0.002219
		self.epc0 = -0.002
		
		m = 1+3600/(self.fpc*ratio)/ratio
		if 'epcc' in kwargs:
			self.epcc = kwargs['ecp']
		else:
			self.epcc = self.epc0*(1+5*(-self.fpcc/self.fpc - 1))
		
		self.Esec = self.fpcc / self.epcc
		self.r_conf = self.Ec/(self.Ec - self.Esec)
		#self.ecu = self.epcc*(1+20*(self.fple/self.fpc))
		self.ecu = - (0.004 + 1.4 * 0.00831 * 68 * 0.09 / 6.899)
		# Debug Variables
		#print(self.__dict__)
		
		self.min_strain = self.ecu
		self.max_strain = 0.002
		
	def stress(self, ec):
		"""Output confined concrete stress for a given strain"""
		self.fail = False
		if ec > 0:
			# Tension Branch
			self.state = 'Cyan'
			#return 0 #TESTING
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
				self.fail = f"Confined Concrete Crushing\nMax Available Strain={round(self.min_strain,5)}\nStrain Experienced={round(ec,5)}"
				return 0


class steel_material:
	"""Model behavior from given steel material properties"""
	def __init__(self, E, fy, fsu, e_sh, e_su, P):
		# Needed steel properties
		self.Es = E
		self.fy = fy
		self.fsu = fsu
		self.e_y = fy/E
		self.e_sh = e_sh
		self.e_su = e_su
		self.P = P
		
		self.min_strain = -self.e_su
		self.max_strain = self.e_su
		
	def stress(self, e):
		"""Output stress in steel from a given strain"""
		self.fail = False
		# TODO - Combine ten/comp, multiply by strain/abs(strain) for sign
		if e < 0:
			#Mirroed Response of tension
			if e > -self.e_y:
				# Elastic
				self.state = 'Yellow'
				return self.Es * e
			elif e > -self.e_sh:
				# Psudo Plastic Before SH
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
				self.fail = f"Steel Fracture\nMax Available Strain={abs(round(self.min_strain,5))}\nStrain Experienced={abs(round(e,5))}"
				return 0
		else:
			### Tension Branch
			if e < self.e_y:
				# Elastic
				self.state = 'Yellow'
				return self.Es * e
			elif e < self.e_sh:
				# Psudo Plastic Before SH
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
				self.fail = f"Steel Fracture\nMax Available Strain={abs(round(self.min_strain,5))}\nStrain Experienced={abs(round(e,5))}"
				return 0

class user_material:
	"""Model behavior from given (strain,stress) points"""
	def __init__(self, points, mirror=False):
		self.strains  = [i[0] for i in points]
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
				sys.exit()
			# if there's a point of 0 strain, don't copy it
			# or maybe there's a np function to remove duplicates?
			if min(self.strains) == 0:
				self.strains  = [-i for i in self.strains[::-1]] + self.strains[1:]
				self.stresses = [-i for i in self.stresses[::-1]] + self.stresses[1:]
			else:
				self.strains  = [-i for i in self.strains[::-1]] + self.strains
				self.stresses = [-i for i in self.stresses[::-1]] + self.stresses
			
		self.min_strain = min(self.strains)
		self.max_strain = max(self.strains)
			
	def stress(self, e):
		"""Output stress from a given strain"""
		self.fail = False
		if (e < min(self.strains)) or (e > max(self.strains)):
			# No more strength - below min or above max strain
			if (self.mirror):
				# stop analysis if user material is mirrored
				self.fail = f"User Material Failure\nMax Available Strain={abs(round(self.min_strain,5))}\nStrain Experienced={abs(round(e,5))}"
			self.state = 'Black'
			return 0
		if abs(e) < self.yield_strain:
			self.state = 6
		else:
			self.state = 7
		# Todo: Want to color with gradient?
		# if self.gradient:
		#	loc = abs(e) / max(max(self.strains),-min(self.strains))
		#	self.state = (loc,0,1-loc)
		return np.interp(e, self.strains, self.stresses)
		
	