# Module Imports
from osfyber.analysis import fiber_model, fiber, Secant_Method, state
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
# Gif Export of disp_mc
import matplotlib.animation as anim

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
        self.state_id = 0
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
    
    def set_load(self, type, **kwargs):
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
        # A = abs(x1*y2 + x2*y3 + x3*y1 - y1*x2 - y2*x3 - y3*x1)/2
        for i,e in enumerate(self.mesh.elements):
            p1 = self.mesh.points[e[0]]
            p2 = self.mesh.points[e[1]]
            p3 = self.mesh.points[e[2]]
            self.mesh_areas[i] = abs(p1[0]*p2[1]+p2[0]*p3[1]+p3[0]*p1[1]-p1[1]*p2[0]-p2[1]*p3[0]-p3[1]*p1[0])/2
        # Assign material ids to elements
        # A bit verbose just to show calcs - might change later (since this only accounts for circle assignements)
        for primative in self.ele_mat_primative:
            R = primative[0]
            prim_c = primative[1]
            mat_id = primative[2]
            for n,c in self.mesh_centroids.items():
                if hypot(c[0]-prim_c[0], c[1]-prim_c[1]) < R:
                    self.ele_mat[n] = mat_id
    
    def gen_fiber_model(self):
        # We can re-mesh our geometry, but our reinforcement is just points
        # So when we add the reinforcement we should also make holes in mesh (TODO)
        max_y = 0
        min_y = 0
        for n in self.mesh_centroids.keys():
            centroid = self.mesh_centroids[n]
            area = self.mesh_areas[n]
            self.analysis_model.fibers[n] = fiber(area, centroid, self.ele_mat[n])
            max_y = max(max_y,centroid[1])
            min_y = max(max_y,centroid[1])
        self.analysis_model.maxy = max_y
        self.analysis_model.miny = min_y
        meh = max(self.mesh_centroids.keys())
        n = 1
        for reinf_group in self.reinforcement:
            for centroid in reinf_group.points:
                # TODO - SHOULD I HAVE PATCHES REGARDLESS OF MAT/REINF?
                # MIGHT BE PROBLEM IF NEED VALUES FROM A SINGLE BAR
                # MAYBE SET MAX VOLUME SEPARATELY?
                self.analysis_model.fibers[n+meh] = fiber(reinf_group.area, centroid, reinf_group.mat_id)
                n += 1
        self.analysis_model.materials = self.materials
        # Generate Inital States
        self.states[0] = state([0]*len(self.mesh.elements)+[8]*(n-1), 0, 0) 
        self.phi_list = [0]
        self.M_list = [0]
    
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
            self.states[step] = self.analysis_model.state
            #self.state_yloc[step] = self.analysis_model.states
            if self.analysis_model.fail:
                print(f"Analysis ended at phi={round(phi,8)}\nFailure: {self.analysis_model.fail}")
                break
            # TODO Not sure If I need - might be good if we don't see a failure condition?
            #if M == 0:
            #    print("M Failure")
            #    break
    
    def color_from_state(self, state):
        """Not sure how I want to store the states for now, so i'll just hold the options in this function"""
        colors = ["Gray","Purple","Pink","Blue","Black","Green","Yellow","Red","White"]
        if type(state) == int:
            return colors[state]
        else:
            return state
    
    def display_material(self, mat_id, loc=None):
        # Plot a specific material with tag mat_id
        material = self.materials[mat_id]
        tmp_strains = material.useful_points
        useful_strians = []
        for point in tmp_strains:
            useful_strians.append(point-1e-6)
            useful_strians.append(point)
            useful_strians.append(point+1e-6)
        lin_strains = np.linspace(material.useful_points[0], material.useful_points[-1], 30)
        strains = np.sort(np.r_[lin_strains, useful_strians[1:-1]])
        stresses = []
        colors = []
        for strain in strains:
            stresses.append(material.stress(strain))
            colors.append(self.color_from_state(material.state))
        stresses = np.array(stresses)
        stresses[ stresses==0 ] = np.nan
        fig,ax = plt.subplots()
        # Setup Plot - Discrete between colors
        i = 0
        x = []
        y = []
        color = colors[0]
        while True:
            new_color = colors[i]
            x.append(strains[i])
            y.append(stresses[i])
            if new_color != color or i == len(strains)-1:
                ax.plot(x, y, c=color)
                color = new_color
                x = [strains[i]]
                y = [stresses[i]]
                if i == len(strains)-1:
                    break
                i -= 1
            i += 1
        # Axis Options
        if loc != None:
            # Plot some location
            ax.plot(loc[0],loc[1],'ro')
        ax.set_title("Material Stress Strain id={}".format(mat_id))
        ax.set_xlabel('Strain (in/in)')
        ax.set_ylabel('Stress (ksi)')
        ax.grid()
    
    def display_materials(self, mat_id=None, loc=None):
        """Plot all loaded materials to verify stress strain curves"""
        if mat_id == None:
            # Plot standard -> all material plots
            for mat_id,material in self.materials.items():
                self.display_material(mat_id, loc)
        else:
            self.display_material(mat_id, loc)
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
                new_patch = CirclePolygon(centroids,reinf_group.radius,8,fc='White',ec='black',zorder=15)
                patches.append(new_patch)
                ax.add_patch(new_patch)
        # Plot nodal points of mesh (useful to auto set axes limts)
        ax.scatter(*zip(*self.mesh.points),c='black',s=4,zorder=5)
        plt.show()

    def mc_gif(self):
        """Testing saving as gif"""
        fig, (ax, ax_mc) = plt.subplots(1, 2, figsize=(self.figsize*2,self.figsize))
        polygons = []
        for e in self.mesh.elements:
            poly_pts = [self.mesh.points[i] for i in e]
            polygons.append(poly_pts)
        patches = []
        # Draw Polygons for mesh and reinformcent
        i = 0
        for p in polygons:
            new_patch = plt.Polygon(p,fc='Grey',ec='Black',zorder=0, picker=.01, gid=str(i))
            patches.append(new_patch)
            ax.add_patch(new_patch)
            i += 1
        if self.reinforcement:
            for reinf_group in self.reinforcement:
                for centroid in reinf_group.points:
                    new_patch = CirclePolygon(centroid,reinf_group.radius,8,fc='White',ec='Black',zorder=15, picker=.01, gid=str(i))
                    patches.append(new_patch)
                    ax.add_patch(new_patch)
                    i += 1
        # Plot nodal points of mesh (useful to auto set axes limts)
        ax.scatter(*zip(*self.mesh.points),c='black',s=4,zorder=5)
        num_steps = len(self.states)-1
        # LEFT MC Plot
        ax_mc.plot(self.phi_list, self.M_list, 'k')
        ax.axis('equal')
        # Location indicator as patch so we can update it easily
        # Ellipse b/c we don't have equal axis ranges for phi/M
        point = Ellipse((0,0),max(self.phi_list)/50,max(self.M_list)/50,fc='red',zorder=10)
        ax_mc.add_patch(point)
        # Scientific notation for phi values
        ax_mc.ticklabel_format(axis='x', style='sci', scilimits=(-2,2), useMathText=True)
        
        # Setup Sliders for interatablility to look through results
        def update(state):
            state = int(state)
            self.state_id = int(state)
            # Update Patch  colors
            state_colors = [self.color_from_state(e) for e in self.states[state].mat_state]
            for patch,color in zip(patches,state_colors):
                patch.set_facecolor(color)
            point.center = (self.phi_list[state],self.M_list[state])
            #fig.canvas.draw_idle()
        ax_slider = plt.axes([0.117, 0.01, 0.79, 0.02], facecolor='White')
        slider_fct = Slider(ax_slider, 'STEP', 0, num_steps, valinit=0, valstep=1, valfmt='%i')
        slider_fct.on_changed(update)
        
        def frame(state):
            # Set State - Updates everything
            slider_fct.set_val(state)
        
        print("Saving animated gif. This will take a while and take several dozen megs of space.\nIt's saving full frames and not optimizing.")
        test = anim.FuncAnimation(fig, frame, frames=range(len(self.states)))
        test.save("Moment_Curvature.gif", fps=8)
        
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
        i = 0
        for p in polygons:
            new_patch = plt.Polygon(p,fc=(0.8,0.8,0.8),ec='black',zorder=0, picker=.01, gid=str(i))
            patches.append(new_patch)
            ax.add_patch(new_patch)
            i += 1
        if self.reinforcement:
            for reinf_group in self.reinforcement:
                for centroid in reinf_group.points:
                    new_patch = CirclePolygon(centroid,reinf_group.radius,8,fc='White',ec='black',zorder=15, picker=.01, gid=str(i))
                    patches.append(new_patch)
                    ax.add_patch(new_patch)
                    i += 1
        # Plot nodal points of mesh (useful to auto set axes limts)
        ax.scatter(*zip(*self.mesh.points),c='black',s=4,zorder=5)
        num_steps = len(self.states)-1
        # LEFT MC Plot
        ax_mc.plot(self.phi_list, self.M_list, 'k')
        ax.axis('equal')
        # Location indicator as patch so we can update it easily
        # Ellipse b/c we don't have equal axis ranges for phi/M
        point = Ellipse((0,0),max(self.phi_list)/50,max(self.M_list)/50,fc='red',zorder=10)
        ax_mc.add_patch(point)
        # Scientific notation for phi values
        ax_mc.ticklabel_format(axis='x', style='sci', scilimits=(-2,2), useMathText=True)
        
        # Setup Sliders for interatablility to look through results
        def update(state):
            state = int(state)
            self.state_id = int(state)
            # Update Patch  colors
            state_colors = [self.color_from_state(e) for e in self.states[state].mat_state]
            for patch,color in zip(patches,state_colors):
                patch.set_facecolor(color)
            point.center = (self.phi_list[state],self.M_list[state])
            fig.canvas.draw_idle()
        ax_slider = plt.axes([0.117, 0.01, 0.79, 0.02], facecolor='White')
        slider_fct = Slider(ax_slider, 'STEP', 0, num_steps, valinit=0, valstep=1, valfmt='%i')
        slider_fct.on_changed(update)
        
        # Setup mesh onclick event to view stress/strain location of each fiber
        def onpick(event):
            #print(event.mouseevent.__dict__)
            #print(event.artist.__dict__)
            #print(event.canvas.__dict__)
            patch_id = int(event.artist._gid)
            strain = self.states[self.state_id].strains[patch_id]
            stress = self.states[self.state_id].stresses[patch_id]
            self.display_materials(self.ele_mat[patch_id], (strain,stress))
            return True
        fig.canvas.mpl_connect('pick_event', onpick)
        
        plt.show()
        
    def display_mc_2x2(self):
        """Plot the moment curvature, and interactive section viewer"""
        # TODO SEPARATE PLOT FOR FORCE/STRESS DISTRIBUTION?
        # TODO PATCH COLLECTION FOR INSTANTANEOUS? COLOR UPDATE?
        fig, axes  = plt.subplots(2, 2, figsize=(self.figsize*1.5,self.figsize*1.5))
        ax             = axes[0,0]
        ax_mc        = axes[0,1]
        ax_strain    = axes[1,0]
        ax_stress    = axes[1,1]
        polygons = []
        for e in self.mesh.elements:
            poly_pts = [self.mesh.points[i] for i in e]
            polygons.append(poly_pts)
        patches = []
        # Draw Polygons for mesh and reinformcent
        for p in polygons:
            new_patch = plt.Polygon(p,fc=(0.8,0.8,0.8),ec='black',zorder=0, picker=.01)
            patches.append(new_patch)
            ax.add_patch(new_patch)
        if self.reinforcement:
            for reinf_group in self.reinforcement:
                for centroid in reinf_group.points:
                    new_patch = CirclePolygon(centroid,reinf_group.radius,8,fc='White',ec='black',zorder=15, picker=.01)
                    patches.append(new_patch)
                    ax.add_patch(new_patch)
        # Plot nodal points of mesh (useful to auto set axes limts)
        ax.scatter(*zip(*self.mesh.points),c='black',s=4,zorder=5)
        num_steps = len(self.states)-1
        # LEFT MC Plot
        ax_mc.plot(self.phi_list, self.M_list, 'k')
        ax.axis('equal')
        # Location indicator as patch so we can update it easily
        # Ellipse b/c we don't have equal axis ranges for phi/M
        point = Ellipse((0,0),max(self.phi_list)/50,max(self.M_list)/50,fc='red',zorder=10)
        ax_mc.add_patch(point)
        # Scientific notation for phi values
        ax_mc.ticklabel_format(axis='x', style='sci', scilimits=(-2,2), useMathText=True)
        
        # Setup Sliders for interatablility to look through results
        def update(phi_step):
            phi_step = int(phi_step)
            self.state_id = int(phi_step)
            # Update Patch  colors
            state_colors = [self.color_from_state(e) for e in self.states[state].mat_state]
            for patch,color in zip(patches,state_colors):
                patch.set_facecolor(color)
            point.center = (self.phi_list[phi_step],self.M_list[phi_step])
            ax_strain.plot([self.states[phi_step].min_strain,self.states[phi_step].max_strain],[0,1],'r')
            fig.canvas.draw_idle()
        ax_slider = plt.axes([0.117, 0.01, 0.79, 0.02], facecolor='white')
        slider_fct = Slider(ax_slider, 'STEP', 0, num_steps, valinit=0, valstep=1, valfmt='%i')
        slider_fct.on_changed(update)
        
        # Debugging strain/stress plots
        ax_strain.plot([0,0],[0,1],'k')
        ax_strain.plot([0,0],[0,1],'r')
        ax_stress.plot([0,1],[0,1])
        
        plt.tight_layout(rect=(0,.05,1,1))
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
            
