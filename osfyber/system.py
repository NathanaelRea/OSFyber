# Module Imports
from osfyber.analysis import fiber_model, fiber, Secant_Method
from osfyber.materials import *
# Math
import numpy as np
from math import hypot
# Meshing
from meshpy.geometry import GeometryBuilder, make_circle, make_box
from meshpy.triangle import MeshInfo, build
# Plotting
import matplotlib.pyplot as plt
from matplotlib.patches import CirclePolygon, Ellipse
from matplotlib.widgets import Slider

class FyberModel:
	"""
	Full Analysis of Fiber Model for a Column
	"""
	def __init__(self, figsize=6):
		# Material Properties
		self.materials = {}
		# Mesh Properties
		self.builder = GeometryBuilder()
		self.reinforcement = []
		self.outer_facets = []
		self.mesh_areas = {}
		self.mesh_centroids = {}
		# Data structure to asign material ids after mesh generation
		self.ele_mat_primative = []
		self.ele_mat = {}
		# Load step Data
		self.load_angle = 0
		self.phi_list = []
		self.M_list = []
		# States of all elements at each load step from analysis
		self.states = {}
		self.state_yloc = {}
		self.analysis_model = fiber_model()
		# Failure check for end of analysis (conf crush/ bar fracture)
		self.fail = False
		# Plot settings
		self.figsize = figsize
	
	def add_material(self, id, type, **kwargs):
		"""Add materials for unconfined concrete, confined concrete, or steel"""
		if type == 'concrete':
			if 'fple' in kwargs:
				self.materials[id] = conf_conc_material(kwargs)
			else:
				self.materials[id] = unconf_conc_material(kwargs)
		elif type == 'steel':
			self.materials[id] = steel_material(kwargs['E'], kwargs['fy'], kwargs['fsu'], 
				kwargs['e_sh'], kwargs['e_su'], kwargs['P'])
		elif type == 'user':
			if 'mirror' in kwargs:
				self.materials[id] = user_material(kwargs['points'],kwargs['mirror'])
			else:
				self.materials[id] = user_material(kwargs['points'])
		else:
			Exit("'{}' is not a material type I recognize".format(type))
	
	def conf_pressure(self, type, **kwargs):
		if type == 'circle':
			return conf_pressure_circ(kwargs['fyh'], kwargs['bar'], kwargs['s'], kwargs['D'])
		else:
			return conf_pressure_rect(kwargs['fyh'], kwargs['bar'], kwargs['s'], kwargs['b'], kwargs['w'], kwargs['nx'], kwargs['ny'])
	
	def add_geometry(self, type, mat_id, **kwargs):
		if type == 'circle':
			self.builder.add_geometry(*make_circle(kwargs['D']/2,kwargs['c']))
			self.ele_mat_primative.append((kwargs['D']/2,kwargs['c'],mat_id))
			if not self.outer_facets:
				# list empty - first points
				self.outer_facets = self.builder.facets
			#self.mesh_points.extend(builder.points)
		elif type == 'rectangle':
			# TODO CHANGE TO builder.add_geometry(*make_rectangle)
			# TODO ADD RECT PRIMIATIVE
			#self.mesh_points.extend(gen_rect_points(c, W, H))
			print("I can only currenly add type='circle' reinforcement")
			sys.exit()
		else:
			print("'{}' is not a geometry type I recognize".format(type))
			sys.exit()

	def add_reinforcement(self, type, mat_id, **kwargs):
		"""Add reinforcement in a circle (for now), confine concrete inside if needed"""
		if type == 'circle':
			if 'conf_id' in kwargs:
				reinf_builder = GeometryBuilder()
				reinf_builder.add_geometry(*make_circle(kwargs['D']/2, kwargs['c'], kwargs['count']))
				r = reinf_props(reinf_builder.points, kwargs['bar'], mat_id)
				self.reinforcement.append(r)
				# Add points for confinement into mesh
				self.builder.add_geometry(*make_circle(kwargs['D']/2,kwargs['c']))
				self.ele_mat_primative.append((kwargs['D']/2,kwargs['c'],kwargs['conf_id']))
			else:
				reinf_builder = GeometryBuilder()
				reinf_builder.add_geometry(*make_circle(kwargs['D']/2, kwargs['c'], kwargs['count']))
				r = reinf_props(reinf_builder.points, kwargs['bar'], mat_id)
				self.reinforcement.append(r)
		else:
			# Rectangular reinforcement
			# not implemented for now
			print("I can only currenly add type='circle' reinforcement")
			sys.exit()
	
	def add_load(self, type, **kwargs):
		"""Add axial load and moment increment/angle"""
		if type == 'Axial':
			self.analysis_model.P = kwargs['P']
		elif type == 'Moment':
			self.analysis_model.M = kwargs['M']
		else:
			print("This is not a load type I recognize")
			sys.exit()
	
	def generate_mesh(self, mesh_size=5):
		"""Setup Mesh Object and Properties
		Go through generation Builder -> MeshInfo -> Mesh"""
		# TODO Look into section-properties on github?
		# TODO Are the values from meshpy.triangle accurate?
		mesh_info = MeshInfo()
		self.builder.set(mesh_info)
		self.mesh = build(mesh_info,max_volume=mesh_size)
		# Calculate Elment Centroids
		# C = [(x1+x2+x3)/3 , (y1+y2+y3)/3]
		for i,e in enumerate(self.mesh.elements):
			p1 = self.mesh.points[e[0]]
			p2 = self.mesh.points[e[1]]
			p3 = self.mesh.points[e[2]]
			self.mesh_centroids[i] = ((p1[0]+p2[0]+p3[0])/3,(p1[1]+p2[1]+p3[1])/3)
		# Calculate Element Areas
		# A = abs(x1*y2+x2*y3+x3*y1-y1*x2-y2*x3-y3*x1)/2
		for i,e in enumerate(self.mesh.elements):
			p1 = self.mesh.points[e[0]]
			p2 = self.mesh.points[e[1]]
			p3 = self.mesh.points[e[2]]
			self.mesh_areas[i] = abs(p1[0]*p2[1]+p2[0]*p3[1]+p3[0]*p1[1]-p1[1]*p2[0]-p2[1]*p3[0]-p3[1]*p1[0])/2
		# Assign material ids to elements
		# set each element to mat id -1
		for n in self.mesh_centroids.keys():
			self.ele_mat[n] = -1
		# A bit verbose just to show calcs - might change later (since this only accounts for circle assignements)
		for primative in self.ele_mat_primative:
			R = primative[0]
			C = primative[1]
			mat_id = primative[2]
			for n,c in self.mesh_centroids.items():
				if hypot(c[0]-C[0], c[1]-C[1]) < R:
					self.ele_mat[n] = mat_id
		# Generate Inital States
		self.states[0] = [0]*len(self.mesh.elements)
		self.phi_list = [0]
		self.M_list = [0]
	
	def gen_fiber_model(self):
		# We can re-mesh our geometry, but our reinforcement is just points
		# So when we add the reinforcement we should also make holes in mesh (TODO)
		for n in self.mesh_centroids.keys():
			centroid = self.mesh_centroids[n]
			area = self.mesh_areas[n]
			self.analysis_model.fibers[n] = fiber(area, centroid, self.ele_mat[n])
		meh = max(self.mesh_centroids.keys())
		n = 1
		for reinf_group in self.reinforcement:
			for centroid in reinf_group.points:
				# TODO - SHOULD I HAVE PATCHES REGARDLESS OF MAT/REINF?
				# MIGHT BE PROBLEM IF NEED VALUES FROM A SINGLE BAR
				# MAYBE SET MAX VOLUME SEPARATELY?
				self.analysis_model.fibers[n+meh] = fiber(reinf_group.area, centroid, reinf_group.mat_id)
				n += 1
		# Save material props into the analysis model
		# TODO - CAN I JUST KEEP THEM IN THE ANALYSIS MODEL FROM THE START?
		self.analysis_model.materials = self.materials
	
	def analyze(self):
		"""Analyze a fiber model with given load step"""
		# setup the model with the current mesh discretization
		self.gen_fiber_model()
		# TODO allow step change?
		delta_phi = 2e-5
		for step in range(1,1000):
			# Setup this step to run equilibrium analysis
			phi = step * delta_phi
			self.analysis_model.phi = phi
			zero_strain = Secant_Method(self.analysis_model.force_balance, self.analysis_model.zero_strain_location)
			self.analysis_model.zero_strain_location = zero_strain
			# Get Moment at current state
			M = self.analysis_model.calc_M()
			# Generate state data and save this load step results
			self.phi_list.append(round(phi,8))
			self.M_list.append(M)
			self.states[step] = self.analysis_model.states
			self.state_yloc[step] = self.analysis_model.states
			if self.analysis_model.fail:
				print(f"Analysis ended at phi={round(phi,8)}\n{self.analysis_model.fail}")
				break
			# TODO Not sure If I need - might be good if we don't see a failure condition?
			#if M == 0:
			#	print("M Failure")
			#	break
	
	def display_mat(self, all=False):
		"""Plot all loaded materials to verify stress strain curves"""
		if all:
			fig,ax = plt.subplots()
			ax.set_xlabel('Strain (in/in)')
			ax.set_ylabel('Stress (ksi)')
			ax.grid()
			for mat_id,material in self.materials.items():
				strains = np.linspace(material.min_strain,material.max_strain,500)
				stress = np.array([material.stress(strain) for strain in strains])
				stress[ stress==0 ] = np.nan
				ax.plot(strains,stress,'k')
		else:
			for mat_id,material in self.materials.items():
				strains = np.linspace(material.min_strain,material.max_strain,500)
				stress = np.array([material.stress(strain) for strain in strains])
				stress[ stress==0 ] = np.nan
				fig,ax = plt.subplots()
				ax.plot(strains,stress,'k')
				ax.set_title("Material Stress Strain id={}".format(mat_id))
				ax.set_xlabel('Strain (in/in)')
				ax.set_ylabel('Stress (ksi)')
				ax.grid()
		plt.show()

	def display_mesh(self):
		"""Show the current geometry and mesh (equivilent to step zero of display_mc)"""
		fig, ax = plt.subplots(figsize=(self.figsize,self.figsize))
		patches = []
		for i,pt_ids in enumerate(self.mesh.elements):
			p = [self.mesh.points[pt_id] for pt_id in pt_ids]
			new_patch = plt.Polygon(p,fc=(0.8,0.8,0.8),ec='black',zorder=0)
			patches.append(new_patch)
			ax.add_patch(new_patch)
		for reinf_group in self.reinforcement:
			for centroids in reinf_group.points:
				new_patch = CirclePolygon(centroids,reinf_group.radius,8,fc=(0.1,0.1,0.1),ec='black',zorder=15)
				patches.append(new_patch)
				ax.add_patch(new_patch)
		# Plot nodal points of mesh (useful to auto set axes limts)
		ax.scatter(*zip(*self.mesh.points),c='black',s=4,zorder=5)
		plt.show()

	def display_mc(self):
		"""Plot the moment curvature, and interactive section viewer"""
		# TODO SEPARATE PLOT FOR FORCE/STRESS DISTRIBUTION?
		# TODO PATCH COLLECTION FOR INSTANTANEOUS? COLOR UPDATE?
		fig, (ax, ax_mc) = plt.subplots(1, 2, figsize=(self.figsize*2,self.figsize))
		polygons = []
		for e in self.mesh.elements:
			poly_pts = [self.mesh.points[i] for i in e]
			polygons.append(poly_pts)
		patches = []
		# Draw Polygons for mesh and reinformcent
		for p in polygons:
			new_patch = plt.Polygon(p,fc=(0.8,0.8,0.8),ec='black',zorder=0)
			patches.append(new_patch)
			ax.add_patch(new_patch)
		if self.reinforcement:
			for reinf_group in self.reinforcement:
				for centroid in reinf_group.points:
					new_patch = CirclePolygon(centroid,reinf_group.radius,8,fc=(0.1,0.1,0.1),ec='black',zorder=15)
					patches.append(new_patch)
					ax.add_patch(new_patch)
		# Plot nodal points of mesh (useful to auto set axes limts)
		ax.scatter(*zip(*self.mesh.points),c='black',s=4,zorder=5)
		## STATE COLOR INFORMATION
		#				 unconf <> y	 conf <> y	   reinf <> y	  usr <> y        Fail
		state_colors = ["Gray","Purple","Pink","Blue","Black","Green","Yellow","Red","White"]
		# TODO ADD 3 COLORS PER MATERIAL
		# Before Yield | Before Max Str | Before Failure
		num_steps = len(self.states)-1
		# Setup Sliders for interactablity
		def update(state):
			state = int(state)
			# Update Concrete Mesh Polygons
			for patch,color in zip(patches,self.states[state]):
				if type(color) == int:
					patch.set_facecolor(state_colors[color])
				else:
					patch.set_facecolor(color)
			point.center = (self.phi_list[state],self.M_list[state])
			fig.canvas.draw_idle()
		ax_slider = plt.axes([0.117, 0.01, 0.79, 0.02], facecolor='white')
		slider_fct = Slider(ax_slider, 'STEP', 0, num_steps, valinit=0, valstep=1)
		slider_fct.on_changed(update)
		# LEFT MC Plot
		ax_mc.plot(self.phi_list, self.M_list, 'k')
		ax.axis('equal')
		# Location indicator as patch so we can update it easily
		# Ellipse b/c we don't have equal axis ranges for phi/M
		point = Ellipse((0,0),max(self.phi_list)/50,max(self.M_list)/50,fc='red',zorder=10)
		ax_mc.add_patch(point)
		# Scientific notation
		ax_mc.set_xticklabels(['{:2.1e}'.format(phi) for phi in self.phi_list], rotation=30)
		plt.show()
		
	def export_results(self, filename=''):
		out = ["phi,M"]
		for phi,M in zip(self.phi_list,self.M_list):
			out.append(f"{phi},{M}")
		if filename == '':
			print("\n".join(out))
		else:
			with open(filename,'w') as f:
				f.write("\n".join(out))
			
