import sys

def Secant_Method(fct, _value, tol=10**-2, _iter=10**-10):
	"""
	Mimimize a given function to zero with the Secant method
	I'm not sure how to find the derivaitve of a black box function,
	And so we can't directly find the derivative, we must first find
	The linear slope between two evaluation points
	fct:	an actual python function with one input/output 
		This is set up in advance to minimize to zero
	start:	the first input
	tol:	The tolerance by witch to quit when the function return
		is < |fct ret|
	iter:	How much to step by in each loop
	"""
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
			return Bisection_Method(fct, tol)
		_value -= cur / slope
	# Secant method should be able to converge within 4 or so steps.
	return Bisection_Method(fct, initial_val, tol)

def Bisection_Method(fct, i_val, tol):
	"""
	Sometimes, the secant/Newton Rhapshon method can deconverge
	This happens when the inial guess is not within the zone of influence
	Bisection is guarentted to always converge, since we start with two
	sides of the function that return pos and neg values
	"""
	# Find a scope - Function values on either side have opposite signs
	scope = 5
	while True:
		# Check if we get +/- on either side of the scope
		a = i_val-scope
		b = i_val+scope
		left = fct(a)
		right = fct(b)
		# check if opposite signs
		if (abs(left)/left != abs(right)/right):
			break
		scope *= 2
	# Try 100 iterations
	nmax = 100
	for n in range(1, nmax): # limit iterations to prevent infinite loop
		c = (a + b)/2 # new midpoint
		if abs(c_val - a_val)/2 <= tol: # solution found
			return c
		a_val = fct(a)
		c_val = fct(c)
		# Check if midpoint and left are same sign
		if c_val/abs(c_val) == a_val/abs(a_val):
			a = c
		else:
			b = c # new interval
	print("Bisectional failed, this is really bad")
	sys.exit()

class state:
	"""State information for a phi step"""
	def __init__(self, mat_state, max_strain, min_strain):
		self.mat_state = mat_state # List of state values for each patch
		self.max_strain = max_strain # Max strain in section
		self.min_strain = min_strain # Min strain in section
		self.strains = [] # List of strain values for each patch
		self.stresses = [] # List of strain values for each patch
		# TODO REMOVE MAX_STRAIN/MIN_STRAIN replace with max(self.strains)

class fiber:
	"""
	A fiber is just a small piece of geometry.
	It's only properties are area, coords, and material.
	"""
	def __init__(self, area, coords, mat_id):
		self.area = area
		self.xy = coords
		self.mat_id = mat_id
		self.fail = False
		# self.distance = None? # for quick calcs

class fiber_model:
	"""
	A method to put all mesh and reinforcement informatino in one place
	That way we can easily calculate the state of the system, and results
	Such as M, and force ballance at a given curvature level
	"""
	def __init__(self):
		"""Analysis Model Class
		Funcitons to find equilibrium and return state data/ M from phi level"""
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
		
	def calc_strain(self, coords, y_intercept):
		# TODO - ADD CURVATURE ANGLE from +y
		# (Currently only does +y direction)
		strain = self.phi * (y_intercept - coords[1])
		return self.phi * (y_intercept - coords[1])
	
	def force_balance(self, y_intercept):
		"""
		Given the intercept (zero_strain_location) and slope (curvature)
		we can find the forces of each fiber and sum them (to be able to balance)
		"""
		self.zero_strain_location = y_intercept
		sum_F = 0
		for i,fiber in self.fibers.items():
			strain = self.calc_strain(fiber.xy, y_intercept)
			stress = self.materials[fiber.mat_id].stress(strain)
			sum_F += fiber.area * stress
		#print(f"{y_intercept}: {sum_F}")
		return sum_F + self.P
	
	def calc_M(self):
		"""Find the internal moment of the centroid to balance the
		External moment and the internal stress from curvature"""
		# TODO: Change this to save state. (Or have another fct that calls this fct)
		# M is just one property of this solved curvature level
		
		# Go through all elemenets (mesh/reinf) and calculate strain from linear centroid, then stress
		# TODO - Find some transformation value c1*x+c2*y -> distance from strain reversal/centroid
		# BUT ONLY CALCULATE THAT ONCE>>>> THEN JUST CALCULAT THE Sin(theta) OR WHATEVER TO FIND THE STEP
		#	WE DO NOT NEED TO RECALCULATE THIS EVERY TIME, BECAUSE WE HAVE A KNOWN STEP SIZE, AND MOMENT ROTATION DOES NOT CHANGE
		
		# Calculate strain at furthest fibers
		max_strain = self.phi * (self.zero_strain_location - self.maxy)
		min_strain = self.phi * (self.zero_strain_location - self.miny)
		
		sum_M = 0
		mat_states = []
		strains = []
		stresses = []
		for i,fiber in self.fibers.items():
			strain = self.calc_strain(fiber.xy, self.zero_strain_location)
			stress = self.materials[fiber.mat_id].stress(strain)
			strains.append(strain)
			stresses.append(stress)
			sum_M += fiber.area * stress * (self.zero_strain_location-fiber.xy[1])
			mat_states.append(self.materials[fiber.mat_id].state)
			fiber.fail = self.materials[fiber.mat_id].fail
			if fiber.fail:
				rnd_loc  = [round(i,3) for i in fiber.xy]
				self.fail = fiber.fail + f"\n\tMat_id={fiber.mat_id}\n\tLocation={rnd_loc}"
		self.state = state(mat_states, max_strain, min_strain)
		self.state.strains = strains
		self.state.stresses = stresses
		return sum_M
