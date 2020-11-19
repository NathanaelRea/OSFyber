# Module Imports
from osfyber.analysis import FiberModel, Fiber, State
from osfyber.materials import *
# Math
import numpy as np
from math import hypot
from scipy.optimize import newton
# Meshing
from meshpy.geometry import GeometryBuilder, make_circle
# TODO from meshpy.geometry import make_box
from meshpy.triangle import MeshInfo, build
# Plotting
import matplotlib.pyplot as plt
from matplotlib.patches import CirclePolygon, Ellipse
from matplotlib.widgets import Slider
# import matplotlib.lines as mlines
# from matplotlib.path import Path
# Gif Export of display_mc
import matplotlib.animation as anim
# Ignore display_mc_2x2() tight layout warning
# import warnings
from sys import exit

"""
TODO better analysis - find better step size based on size?
TODO better analysis - find maximum tension/compression and throw error if static force too high
"""


class FyberModel:
    """
    Full Fiber Model
    Can store all material and geometric properties
    Can run pushover of model
    Can display results in a variety of ways
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
        self.mesh = None
        # Data structure to assign material ids after mesh generation
        self.ele_mat_primitive = []
        self.ele_mat = {}
        # Load step Data
        self.load_angle = 0
        self.phi_list = []
        self.M_list = []
        # States of all elements at each load step from analysis
        self.state_id = 0
        self.states = {}
        self.state_y_loc = {}
        self.analysis_model = FiberModel()
        # Failure check for end of analysis (conf crush/ bar fracture)
        self.fail = False
        # Plot settings
        self.figsize = figsize

    @staticmethod
    def conf_pressure(shape, **kwargs):
        """Technically static method, but I just want one import"""
        if shape == 'circle':
            return conf_pressure_circle(kwargs['fyh'], kwargs['bar'], kwargs['s'], kwargs['D'])
        elif shape == 'rect':
            return conf_pressure_rect(kwargs['fyh'], kwargs['bar'], kwargs['s'], kwargs['b'],
                                      kwargs['w'], kwargs['nx'], kwargs['ny'])
        else:
            exit("I can only find the confining pressure of 'circle' or 'rect'")

    @staticmethod
    def color_from_state(state):
        """Not sure how I want to store the states for now, so i'll just hold the options in this function"""
        colors = ["Gray", "Purple", "Pink", "Blue", "Black", "Green", "Yellow", "Red", "White"]
        if type(state) == int:
            return colors[state]
        else:
            return state

    def add_material(self, mat_num, mat, **kwargs):
        """Add materials for unconfined concrete, confined concrete, or steel"""
        if mat == 'concrete':
            if 'fple' in kwargs:
                self.materials[mat_num] = ConfConcMat(kwargs)
            else:
                self.materials[mat_num] = UnconfConcMat(kwargs)
        elif mat == 'steel':
            self.materials[mat_num] = SteelMat(kwargs['E'], kwargs['fy'], kwargs['fsu'],
                                               kwargs['e_sh'], kwargs['e_su'], kwargs['P'])
        elif mat == 'user':
            if 'mirror' in kwargs:
                self.materials[mat_num] = UserMat(kwargs['points'], kwargs['mirror'])
            else:
                self.materials[mat_num] = UserMat(kwargs['points'])
        else:
            exit("'{}' is not a material type I recognize".format(mat))

    def add_geometry(self, shape, mat_id, **kwargs):
        if shape == 'circle':
            self.builder.add_geometry(*make_circle(kwargs['D'] / 2, kwargs['c']))
            self.ele_mat_primitive.append((kwargs['D'] / 2, kwargs['c'], mat_id))
            if not self.outer_facets:
                # list empty - first points
                self.outer_facets = self.builder.facets
            # self.mesh_points.extend(builder.points)
        elif shape == 'rectangle':
            # TODO CHANGE TO builder.add_geometry(*make_rectangle)
            # TODO add rect primitive
            # self.mesh_points.extend(gen_rect_points(c, W, H))
            print("I can only currently add type='circle' reinforcement")
            exit()
        else:
            print("'{}' is not a geometry type I recognize".format(shape))
            exit()

    def add_reinforcement(self, shape, mat_id, **kwargs):
        """Add reinforcement in a circle (for now), confine concrete inside if needed"""
        if shape == 'circle':
            if 'conf_id' in kwargs:
                reinforcement_builder = GeometryBuilder()
                reinforcement_builder.add_geometry(*make_circle(kwargs['D'] / 2, kwargs['c'], kwargs['count']))
                r = ReinforcementProperties(reinforcement_builder.points, kwargs['bar'], mat_id)
                self.reinforcement.append(r)
                # Add points for confinement into mesh
                self.builder.add_geometry(*make_circle(kwargs['D'] / 2, kwargs['c']))
                self.ele_mat_primitive.append((kwargs['D'] / 2, kwargs['c'], kwargs['conf_id']))
            else:
                reinforcement_builder = GeometryBuilder()
                reinforcement_builder.add_geometry(*make_circle(kwargs['D'] / 2, kwargs['c'], kwargs['count']))
                r = ReinforcementProperties(reinforcement_builder.points, kwargs['bar'], mat_id)
                self.reinforcement.append(r)
        else:
            # Rectangular reinforcement
            # not implemented for now
            print("I can only currently add type='circle' reinforcement")
            exit()

    def set_load(self, load_type, **kwargs):
        """Add axial load and moment increment/angle"""
        if load_type == 'Axial':
            self.analysis_model.P = kwargs['P']
        elif load_type == 'Moment':
            self.analysis_model.M = kwargs['M']
        else:
            print("This is not a load type I recognize")
            exit()

    def generate_mesh(self, mesh_size=5):
        """Setup Mesh Object and Properties
        Go through generation Builder -> MeshInfo -> Mesh"""
        # TODO Look into section-properties on github?
        # TODO Are the values from meshpy.triangle accurate?
        mesh_info = MeshInfo()
        self.builder.set(mesh_info)
        self.mesh = build(mesh_info, max_volume=mesh_size)
        # Calculate Element Centroids
        # C = [(x1+x2+x3)/3 , (y1+y2+y3)/3]
        for i, e in enumerate(self.mesh.elements):
            p1 = self.mesh.points[e[0]]
            p2 = self.mesh.points[e[1]]
            p3 = self.mesh.points[e[2]]
            self.mesh_centroids[i] = ((p1[0] + p2[0] + p3[0]) / 3, (p1[1] + p2[1] + p3[1]) / 3)
        # Calculate Element Areas
        # A = abs(x1*y2 + x2*y3 + x3*y1 - y1*x2 - y2*x3 - y3*x1)/2
        for i, e in enumerate(self.mesh.elements):
            p1 = self.mesh.points[e[0]]
            p2 = self.mesh.points[e[1]]
            p3 = self.mesh.points[e[2]]
            self.mesh_areas[i] = abs(
                p1[0] * p2[1] + p2[0] * p3[1] + p3[0] * p1[1] - p1[1] * p2[0] - p2[1] * p3[0] - p3[1] * p1[0]) / 2
        # Assign material ids to elements
        # A bit verbose just to show calculations - might change later
        # this only accounts for circle assignments
        for primitive in self.ele_mat_primitive:
            r = primitive[0]
            prim_c = primitive[1]
            mat_id = primitive[2]
            for n, c in self.mesh_centroids.items():
                if hypot(c[0] - prim_c[0], c[1] - prim_c[1]) < r:
                    self.ele_mat[n] = mat_id

    def gen_fiber_model(self):
        # We can re-mesh our geometry, but our reinforcement is just points
        # So when we add the reinforcement we should also make holes in mesh (TODO)
        max_y = 0
        min_y = 0
        for n in self.mesh_centroids.keys():
            centroid = self.mesh_centroids[n]
            area = self.mesh_areas[n]
            self.analysis_model.fibers[n] = Fiber(area, centroid, self.ele_mat[n])
            max_y = max(max_y, centroid[1])
            min_y = min(min_y, centroid[1])
        self.analysis_model.maxy = max_y
        self.analysis_model.miny = min_y
        meh = max(self.mesh_centroids.keys())
        n = 1
        for reinforce_group in self.reinforcement:
            for centroid in reinforce_group.points:
                # TODO - SHOULD I HAVE PATCHES REGARDLESS OF MAT/reinforcement?
                # MIGHT BE PROBLEM IF NEED VALUES FROM A SINGLE BAR
                # MAYBE SET MAX VOLUME SEPARATELY?
                self.analysis_model.fibers[n + meh] = Fiber(reinforce_group.area, centroid, reinforce_group.mat_id)
                n += 1
        self.analysis_model.materials = self.materials
        # Generate Initial States
        self.states[0] = State([0] * len(self.mesh.elements) + [8] * (n - 1), 0, 0, 0)
        self.phi_list = [0]
        self.M_list = [0]

    def analyze(self):
        """Analyze a fiber model with given load step"""
        # setup the model with the current mesh discretization
        self.gen_fiber_model()
        # TODO allow step change?
        # delta_phi = 5e-6
        delta_phi = 1e-5
        for step in range(1, 1000):
            # Setup this step to run equilibrium analysis
            phi = step * delta_phi
            self.analysis_model.phi = phi
            # I'm not sure if scipy.optimize.newton can fail here, might need to catch exception with bisection
            zero_strain = newton(self.analysis_model.force_balance, self.analysis_model.zero_strain_location)
            self.analysis_model.zero_strain_location = zero_strain
            # Get Moment at current state
            moment = self.analysis_model.calc_moment()
            # Generate state data and save this load step results
            self.phi_list.append(round(phi, 8))
            self.M_list.append(moment)
            self.states[step] = self.analysis_model.state
            # self.state_y_loc[step] = self.analysis_model.states
            if self.analysis_model.fail:
                print(f"Analysis ended at phi={round(phi, 8)}\nFailure: {self.analysis_model.fail}")
                break
            # Don't necessarily want to break for tension failure - Say, if entire section is steel
            # if M < 1:
            #   print("Capacity of entire section failed")
            #   break

    def display_material(self, mat_id, loc=None):
        # Plot a specific material with tag mat_id
        material = self.materials[mat_id]
        tmp_strains = material.useful_points
        useful_strains = []
        for point in tmp_strains:
            useful_strains.append(point - 1e-6)
            useful_strains.append(point)
            useful_strains.append(point + 1e-6)
        lin_strains = np.linspace(material.useful_points[0], material.useful_points[-1], 30)
        strains = np.sort(np.r_[lin_strains, useful_strains[1:-1]])
        stresses = []
        colors = []
        for strain in strains:
            stresses.append(material.stress(strain))
            colors.append(self.color_from_state(material.state))
        stresses = np.array(stresses)
        stresses[stresses == 0] = np.nan
        fig, ax = plt.subplots()
        # Setup Plot - Discrete between colors
        i = 0
        x = []
        y = []
        color = colors[0]
        while True:
            new_color = colors[i]
            x.append(strains[i])
            y.append(stresses[i])
            if new_color != color or i == len(strains) - 1:
                ax.plot(x, y, c=color)
                color = new_color
                x = [strains[i]]
                y = [stresses[i]]
                if i == len(strains) - 1:
                    break
                i -= 1
            i += 1
        # Axis Options
        if loc is not None:
            # Plot some location
            ax.plot(loc[0], loc[1], 'ro')
        ax.set_title("Material Stress Strain id={}".format(mat_id))
        ax.set_xlabel('Strain (in/in)')
        ax.set_ylabel('Stress (ksi)')
        ax.grid()

    def display_materials(self, mat_id=None, loc=None):
        """Plot all loaded materials to verify stress strain curves"""
        if mat_id is None:
            # Plot standard -> all material plots
            for mat_id, material in self.materials.items():
                self.display_material(mat_id, loc)
        else:
            self.display_material(mat_id, loc)
        plt.show()

    def display_mesh(self):
        """Show the current geometry and mesh (equivalent to step zero of display_mc)"""
        fig, ax = plt.subplots(figsize=(self.figsize, self.figsize))
        patches = []
        for i, pt_ids in enumerate(self.mesh.elements):
            p = [self.mesh.points[pt_id] for pt_id in pt_ids]
            new_patch = plt.Polygon(p, fc="Gray", ec='black', zorder=0)
            patches.append(new_patch)
            ax.add_patch(new_patch)
        for reinforce_group in self.reinforcement:
            for centroids in reinforce_group.points:
                new_patch = CirclePolygon(centroids, reinforce_group.radius, 8, fc='White', ec='black', zorder=15)
                patches.append(new_patch)
                ax.add_patch(new_patch)
        # Plot nodal points of mesh (useful to auto set axes limits)
        ax.scatter(*zip(*self.mesh.points), c='black', s=4, zorder=5)
        plt.show()

    def mc_gif(self):
        """Testing saving as gif"""
        fig, (ax, ax_mc) = plt.subplots(1, 2, figsize=(self.figsize * 2, self.figsize))
        polygons = []
        for element in self.mesh.elements:
            poly_pts = [self.mesh.points[i] for i in element]
            polygons.append(poly_pts)
        patches = []
        # Draw Polygons for mesh and reinforcement
        i = 0
        for p in polygons:
            new_patch = plt.Polygon(p, fc='Grey', ec='Black', zorder=0, picker=.01, gid=str(i))
            patches.append(new_patch)
            ax.add_patch(new_patch)
            i += 1
        if self.reinforcement:
            for reinforcement_group in self.reinforcement:
                for centroid in reinforcement_group.points:
                    new_patch = CirclePolygon(centroid, reinforcement_group.radius, 8, fc='White', ec='Black',
                                              zorder=15, picker=.01, gid=str(i))
                    patches.append(new_patch)
                    ax.add_patch(new_patch)
                    i += 1
        # Plot nodal points of mesh (useful to auto set axes limits)
        ax.scatter(*zip(*self.mesh.points), c='black', s=4, zorder=5)
        num_steps = len(self.states) - 1
        # LEFT MC Plot
        ax_mc.plot(self.phi_list, self.M_list, 'k')
        ax.axis('equal')
        # Location indicator as patch so we can update it easily
        # Ellipse b/c we don't have equal axis ranges for phi/M
        point = Ellipse((0, 0), max(self.phi_list) / 50, max(self.M_list) / 50, fc='red', zorder=10)
        ax_mc.add_patch(point)
        # Scientific notation for phi values
        ax_mc.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2), useMathText=True)

        # Setup Sliders for intractability to look through results
        def update(state):
            state = int(state)
            self.state_id = int(state)
            # Update Patch  colors
            state_colors = [self.color_from_state(state) for state in self.states[state].mat_state]
            for patch, color in zip(patches, state_colors):
                patch.set_facecolor(color)
            point.center = (self.phi_list[state], self.M_list[state])
            # fig.canvas.draw_idle()

        ax_slider = plt.axes([0.117, 0.01, 0.79, 0.02], facecolor='White')
        slider_fct = Slider(ax_slider, 'STEP', 0, num_steps, valinit=0, valstep=1, valfmt='%i')
        slider_fct.on_changed(update)

        def frame(state):
            # Set State - Updates everything
            slider_fct.set_val(state)

        print("Saving animated gif. This will take a while and take several dozen megs of space.\n"
              "It's saving full frames and not optimizing.")
        test = anim.FuncAnimation(fig, frame, frames=len(self.states))
        test.save("Moment_Curvature.gif", fps=16)

    def display_mc(self):
        """Plot the moment curvature, and interactive section viewer"""
        # TODO SEPARATE PLOT FOR FORCE/STRESS DISTRIBUTION?
        # TODO PATCH COLLECTION FOR INSTANTANEOUS? COLOR UPDATE?
        fig, (ax, ax_mc) = plt.subplots(1, 2, figsize=(self.figsize * 2, self.figsize))
        polygons = []
        for element in self.mesh.elements:
            poly_pts = [self.mesh.points[i] for i in element]
            polygons.append(poly_pts)
        patches = []
        # Draw Polygons for mesh and reinforcement
        i = 0
        for p in polygons:
            new_patch = plt.Polygon(p, fc="Gray", ec='black', zorder=0, picker=.01, gid=str(i))
            patches.append(new_patch)
            ax.add_patch(new_patch)
            i += 1
        if self.reinforcement:
            for reinforce_group in self.reinforcement:
                for centroid in reinforce_group.points:
                    new_patch = CirclePolygon(centroid, reinforce_group.radius, 8, fc='White', ec='black',
                                              zorder=15, picker=.01, gid=str(i))
                    patches.append(new_patch)
                    ax.add_patch(new_patch)
                    i += 1
        # Plot nodal points of mesh (useful to auto set axes limits)
        ax.scatter(*zip(*self.mesh.points), c='black', s=4, zorder=5)
        num_steps = len(self.states) - 1
        # LEFT MC Plot
        ax_mc.plot(self.phi_list, self.M_list, 'k')
        ax.axis('equal')
        # Location indicator as patch so we can update it easily
        # Ellipse b/c we don't have equal axis ranges for phi/M
        point = Ellipse((0, 0), max(self.phi_list) / 50, max(self.M_list) / 50, fc='red', zorder=10)
        ax_mc.add_patch(point)
        # Scientific notation for phi values
        ax_mc.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2), useMathText=True)

        # Setup Sliders for intractability to look through results
        def update(phi_step):
            phi_step = int(phi_step)
            self.state_id = int(phi_step)
            # Update Patch  colors
            state_colors = [self.color_from_state(state) for state in self.states[phi_step].mat_state]
            for patch, color in zip(patches, state_colors):
                patch.set_facecolor(color)
            point.center = (self.phi_list[phi_step], self.M_list[phi_step])
            fig.canvas.draw_idle()

        ax_slider = plt.axes([0.117, 0.01, 0.79, 0.02], facecolor='White')
        slider_fct = Slider(ax_slider, 'STEP', 0, num_steps, valinit=0, valstep=1, valfmt='%i')
        slider_fct.on_changed(update)

        # Setup mesh onclick event to view stress/strain location of each fiber
        def onclick(event):
            # print(event.mouseevent.__dict__)
            # print(event.artist.__dict__)
            # print(event.canvas.__dict__)
            patch_id = int(event.artist._gid)
            strain = self.states[self.state_id].strains[patch_id]
            stress = self.states[self.state_id].stresses[patch_id]
            self.display_materials(self.ele_mat[patch_id], (strain, stress))
            return True

        fig.canvas.mpl_connect('pick_event', onclick)

        plt.show()

    def display_mc_2x2(self):
        """Plot the moment curvature, and interactive section viewer"""
        # TODO SEPARATE PLOT FOR FORCE/STRESS DISTRIBUTION?
        # TODO PATCH COLLECTION FOR INSTANTANEOUS? COLOR UPDATE?
        fig, axes = plt.subplots(2, 2, figsize=(self.figsize * 1.5, self.figsize * 1.5))
        ax = axes[0, 0]
        ax_mc = axes[0, 1]
        ax_strain = axes[1, 0]
        ax_stress = axes[1, 1]
        polygons = []
        for element in self.mesh.elements:
            poly_pts = [self.mesh.points[i] for i in element]
            polygons.append(poly_pts)
        patches = []
        # Draw Polygons for mesh and reinforcement
        i = 0
        for p in polygons:
            new_patch = plt.Polygon(p, fc="Gray", ec='black', zorder=0, picker=.01, gid=str(i))
            patches.append(new_patch)
            ax.add_patch(new_patch)
            i += 1
        if self.reinforcement:
            for reinforce_group in self.reinforcement:
                for centroid in reinforce_group.points:
                    new_patch = CirclePolygon(centroid, reinforce_group.radius, 8, fc='White', ec='black',
                                              zorder=15, picker=.01, gid=str(i))
                    patches.append(new_patch)
                    ax.add_patch(new_patch)
                    i += 1
        # Plot nodal points of mesh (useful to auto set axes limits)
        ax.scatter(*zip(*self.mesh.points), c='black', s=4, zorder=5)
        num_steps = len(self.states) - 1
        # LEFT MC Plot
        ax_mc.plot(self.phi_list, self.M_list, 'k')
        ax.axis('equal')
        # Location indicator as patch so we can update it easily
        # Ellipse b/c we don't have equal axis ranges for phi/M
        point = Ellipse((0, 0), max(self.phi_list) / 50, max(self.M_list) / 50, fc='red', zorder=10)
        ax_mc.add_patch(point)
        # Scientific notation for phi values
        ax_mc.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2), useMathText=True)

        line = line2 = line3 = None

        # Setup Sliders for intractability to look through results
        def update(phi_step):
            phi_step = int(phi_step)
            self.state_id = int(phi_step)
            # Update Patch  colors
            state_colors = [self.color_from_state(state) for state in self.states[phi_step].mat_state]
            for patch, color in zip(patches, state_colors):
                patch.set_facecolor(color)
            point.center = (self.phi_list[phi_step], self.M_list[phi_step])

            x = [self.states[phi_step].min_strain, self.states[phi_step].max_strain]
            y = [self.analysis_model.miny, self.analysis_model.maxy]

            line.set_data(x, y)
            line2.set_data([-1, 1], [self.states[phi_step].y_loc, self.states[phi_step].y_loc])
            line3.set_data([self.analysis_model.miny * 1.1, self.analysis_model.maxy * 1.1],
                           [self.states[phi_step].y_loc, self.states[phi_step].y_loc])
            fig.canvas.draw_idle()

        # Setup mesh onclick event to view stress/strain location of each fiber
        def onclick(event):
            # print(event.mouseevent.__dict__)
            # print(event.artist.__dict__)
            # print(event.canvas.__dict__)
            patch_id = int(event.artist._gid)
            strain = self.states[self.state_id].strains[patch_id]
            stress = self.states[self.state_id].stresses[patch_id]
            self.display_materials(self.ele_mat[patch_id], (strain, stress))
            return True

        fig.canvas.mpl_connect('pick_event', onclick)

        ax_slider = plt.axes([0.117, 0.01, 0.79, 0.02], facecolor='white')
        slider_fct = Slider(ax_slider, 'STEP', 0, num_steps, valinit=0, valstep=1, valfmt='%i')
        slider_fct.on_changed(update)

        # Strain plot
        ax_strain.plot([0, 0], [self.analysis_model.miny, self.analysis_model.maxy], 'k')
        ax_strain.set_xlim([self.states[len(self.states) - 1].min_strain * 1.05,
                            self.states[len(self.states) - 1].max_strain * 1.05])

        line, = ax_strain.plot([0, 0], [self.analysis_model.miny, self.analysis_model.maxy], color='Red')
        phi_step_init = 0
        line2, = ax_strain.plot([-1, 1], [self.states[phi_step_init].y_loc, self.states[phi_step_init].y_loc],
                                color='Blue')
        line3, = ax.plot([self.analysis_model.miny * 1.1, self.analysis_model.maxy * 1.1],
                         [self.states[phi_step_init].y_loc, self.states[phi_step_init].y_loc], color='Red')

        # Stress Plot
        ax_stress.plot([0, 0], [self.analysis_model.miny, self.analysis_model.maxy], 'k')
        ax_stress.text(0, 0, "This Plot Was\nIntentionally Left Blank\n(For Now)")

        plt.tight_layout()
        # plt.tight_layout(rect=[0, .05, 1, 1])
        plt.show()

    def export_results(self, filename=''):
        out = ["phi,M"]
        for phi, M in zip(self.phi_list, self.M_list):
            out.append(f"{phi},{M}")
        if filename == '':
            print("\n".join(out))
        else:
            with open(filename, 'w') as f:
                f.write("\n".join(out))
