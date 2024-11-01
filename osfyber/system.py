from dataclasses import dataclass
from typing import Any, Optional

from matplotlib.backend_bases import Event, PickEvent
from matplotlib.figure import Figure
from osfyber.analysis import Fiber, FiberModel, State
from osfyber.materials import (
    Color,
    ConfConcMat,
    ConfinedConcreteProps,
    Material,
    ReinforcementProperties,
    SteelMat,
    SteelProps,
    UnconfConcMat,
    UnconfinedConcreteProps,
    UserMat,
    UserMaterialProps,
    conf_pressure_circle,
    conf_pressure_rect,
)


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
from matplotlib.patches import CirclePolygon, Ellipse, Patch, RegularPolygon
from matplotlib.widgets import Slider

# import matplotlib.lines as mlines
# from matplotlib.path import Path
# Gif Export of display_mc
import matplotlib.animation as anim

# Ignore display_mc_2x2() tight layout warning
# import warnings
from sys import exit


@dataclass
class CircleGeometryProps:
    D: float
    c: tuple[float, float]


@dataclass
class ReinforcementProps:
    bar: int
    count: int
    conf_id: Optional[int] = None


# @dataclass
# class RectGeometryProps:
#     mat_id: int
#     W: float
#     H: float
#     c: tuple[float, float]


@dataclass
class LoadProps:
    P: Optional[float] = None


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

    def __init__(self, figsize: int = 6):
        # Material Properties
        self.materials: dict[int, Material] = {}
        # Mesh Properties
        self.builder = GeometryBuilder()
        self.reinforcement: list[ReinforcementProperties] = []
        self.outer_facets: list[tuple[float, float]] = []
        self.mesh_areas: dict[int, float] = {}
        self.mesh_centroids: dict[int, tuple[float, float]] = {}
        self.mesh: MeshInfo = None
        # Data structure to assign material ids after mesh generation
        self.ele_mat_primitive: list[tuple[float, tuple[float, float], int]] = []
        self.ele_mat: dict[int, int] = {}
        # Load step Data
        self.phi_list: list[float] = []
        self.M_list: list[float] = []
        # States of all elements at each load step from analysis
        self.state_id = 0
        self.states: dict[int, State] = {}
        self.analysis_model = FiberModel()
        # Failure check for end of analysis (conf crush/ bar fracture)
        self.fail = False
        # Plot settings
        self.figsize = figsize

    def add_material_concrete(
        self, mat_num: int, props: UnconfinedConcreteProps
    ) -> None:
        self.materials[mat_num] = UnconfConcMat(props)

    def add_material_confined_concrete(
        self, mat_num: int, props: ConfinedConcreteProps
    ) -> None:
        self.materials[mat_num] = ConfConcMat(props)

    def add_material_steel(self, mat_num: int, props: SteelProps) -> None:
        self.materials[mat_num] = SteelMat(props)

    def add_material_user(self, mat_num: int, props: UserMaterialProps) -> None:
        self.materials[mat_num] = UserMat(props)

    def add_geometry_circle(self, mat_id: int, props: CircleGeometryProps) -> None:
        self.builder.add_geometry(*make_circle(props.D / 2, props.c))
        self.ele_mat_primitive.append((props.D / 2, props.c, mat_id))
        if not self.outer_facets:
            # list empty - first points
            self.outer_facets = self.builder.facets
        # self.mesh_points.extend(builder.points)

    def add_reinforcement_circle(
        self, mat_id: int, props: CircleGeometryProps, props2: ReinforcementProps
    ) -> None:
        """Add reinforcement in a circle (for now), confine concrete inside if needed"""
        if props2.conf_id:
            reinforcement_builder = GeometryBuilder()
            reinforcement_builder.add_geometry(
                *make_circle(props.D / 2, props.c, props2.count)
            )
            r = ReinforcementProperties(
                reinforcement_builder.points, props2.bar, mat_id
            )
            self.reinforcement.append(r)
            # Add points for confinement into mesh
            self.builder.add_geometry(*make_circle(props.D / 2, props.c))
            self.ele_mat_primitive.append((props.D / 2, props.c, props2.conf_id))
        else:
            reinforcement_builder = GeometryBuilder()
            reinforcement_builder.add_geometry(
                *make_circle(props.D / 2, props.c, props2.count)
            )
            r = ReinforcementProperties(
                reinforcement_builder.points, props2.bar, mat_id
            )
            self.reinforcement.append(r)

    def set_load(self, props: LoadProps) -> None:
        """Add axial load"""
        self.analysis_model.P = props.P if props.P else 0.0

    def generate_mesh(self, mesh_size: int = 5) -> None:
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
            self.mesh_centroids[i] = (
                (p1[0] + p2[0] + p3[0]) / 3,
                (p1[1] + p2[1] + p3[1]) / 3,
            )
        # Calculate Element Areas
        # A = abs(x1*y2 + x2*y3 + x3*y1 - y1*x2 - y2*x3 - y3*x1)/2
        for i, e in enumerate(self.mesh.elements):
            p1 = self.mesh.points[e[0]]
            p2 = self.mesh.points[e[1]]
            p3 = self.mesh.points[e[2]]
            self.mesh_areas[i] = (
                abs(
                    p1[0] * p2[1]
                    + p2[0] * p3[1]
                    + p3[0] * p1[1]
                    - p1[1] * p2[0]
                    - p2[1] * p3[0]
                    - p3[1] * p1[0]
                )
                / 2
            )
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

    def gen_fiber_model(self) -> None:
        # We can re-mesh our geometry, but our reinforcement is just points
        # So when we add the reinforcement we should also make holes in mesh (TODO)
        max_y = 0.0
        min_y = 0.0
        for n in self.mesh_centroids.keys():
            centroid = self.mesh_centroids[n]
            area = self.mesh_areas[n]
            self.analysis_model.fibers[n] = Fiber(area, centroid, self.ele_mat[n])
            max_y = max(max_y, centroid[1])
            min_y = min(min_y, centroid[1])
        self.analysis_model.max_y = max_y
        self.analysis_model.min_y = min_y
        meh = max(self.mesh_centroids.keys())
        n = 1
        for reinforce_group in self.reinforcement:
            for centroid in reinforce_group.points:
                # TODO - SHOULD I HAVE PATCHES REGARDLESS OF MAT/reinforcement?
                # MIGHT BE PROBLEM IF NEED VALUES FROM A SINGLE BAR
                # MAYBE SET MAX VOLUME SEPARATELY?
                self.analysis_model.fibers[n + meh] = Fiber(
                    reinforce_group.area, centroid, reinforce_group.mat_id
                )
                n += 1
        self.analysis_model.materials = self.materials
        # Generate Initial States
        self.states[0] = State(
            [Color.Gray] * len(self.mesh.elements) + [Color.White] * (n - 1), 0, 0, 0
        )
        self.phi_list = [0]
        self.M_list = [0]

    def analyze(self) -> None:
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
            zero_strain = newton(
                self.analysis_model.force_balance,
                self.analysis_model.zero_strain_location,
            )
            self.analysis_model.zero_strain_location = zero_strain
            # Get Moment at current state
            moment = self.analysis_model.calc_moment()
            # Generate state data and save this load step results
            self.phi_list.append(round(phi, 8))
            self.M_list.append(moment)
            # FIXME we initialize in gen_fiber_model, probably don't need to? maybe just append?
            if self.analysis_model.state:
                self.states[step] = self.analysis_model.state
            if self.analysis_model.fail:
                print(
                    f"Analysis ended at phi={round(phi, 8)}\nFailure: {self.analysis_model.fail}"
                )
                break
            # Don't necessarily want to break for tension failure - Say, if entire section is steel
            # if M < 1:
            #   print("Capacity of entire section failed")
            #   break

    def display_material(
        self, mat_id: int, loc: Optional[tuple[float, float]] = None
    ) -> None:
        # Plot a specific material with tag mat_id
        material = self.materials[mat_id]
        tmp_strains = list(material.useful_points)
        useful_strains = []
        for point in tmp_strains:
            useful_strains.append(point - 1e-6)
            useful_strains.append(point)
            useful_strains.append(point + 1e-6)
        lin_strains = np.linspace(
            material.useful_points[0], material.useful_points[-1], 30
        )
        strains = np.sort(np.r_[lin_strains, useful_strains[1:-1]])
        out = [material.stress_state(strain) for strain in strains]
        stresses = np.array([s for s, _ in out])
        colors = [c.name for _, c in out]
        stresses[stresses == 0] = np.nan
        _, ax = plt.subplots()
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
            ax.plot(loc[0], loc[1], "ro")
        ax.set_title("Material Stress Strain id={}".format(mat_id))
        ax.set_xlabel("Strain (in/in)")
        ax.set_ylabel("Stress (ksi)")
        ax.grid()
        plt.show()

    def display_materials(self) -> None:
        for mat_id in self.materials.keys():
            self.display_material(mat_id)
        plt.show()

    def display_mesh(self) -> None:
        """Show the current geometry and mesh (equivalent to step zero of display_mc)"""
        fig, ax = plt.subplots(figsize=(self.figsize, self.figsize))
        patches: list[Patch] = []
        for i, pt_ids in enumerate(self.mesh.elements):
            p = [self.mesh.points[pt_id] for pt_id in pt_ids]
            patch = plt.Polygon(p, fc="Gray", ec="black", zorder=0)
            patches.append(patch)
            ax.add_patch(patch)
        for reinforce_group in self.reinforcement:
            for centroids in reinforce_group.points:
                new_patch = CirclePolygon(centroids, reinforce_group.radius)
                new_patch.set_facecolor("White")
                new_patch.set_edgecolor("Black")
                new_patch.set_zorder(15)
                patches.append(new_patch)
                ax.add_patch(new_patch)
        # Plot nodal points of mesh (useful to auto set axes limits)
        x, y = zip(*self.mesh.points)
        ax.scatter(x, y, c="black", s=4, zorder=5)
        plt.show()

    def color_from_state(self, state: Color) -> str:
        # FIXME shouldn't be needed now with enum
        return state.name

    def mc_gif(self) -> None:
        """Testing saving as gif"""
        fig, (ax, ax_mc) = plt.subplots(1, 2, figsize=(self.figsize * 2, self.figsize))
        polygons = []
        for element in self.mesh.elements:
            poly_pts = [self.mesh.points[i] for i in element]
            polygons.append(poly_pts)
        patches: list[Patch] = []
        # Draw Polygons for mesh and reinforcement
        i = 0
        for p in polygons:
            patch = plt.Polygon(
                p, fc="Grey", ec="Black", zorder=0, picker=0.01, gid=str(i)
            )
            patches.append(patch)
            ax.add_patch(patch)
            i += 1
        if self.reinforcement:
            for reinforcement_group in self.reinforcement:
                for centroid in reinforcement_group.points:
                    new_patch = CirclePolygon(centroid, reinforcement_group.radius)
                    new_patch.set_facecolor("White")
                    new_patch.set_edgecolor("Black")
                    new_patch.set_picker(0.01)
                    new_patch.set_gid(str(i))
                    new_patch.set_zorder(15)
                    patches.append(new_patch)
                    ax.add_patch(new_patch)
                    i += 1
        # Plot nodal points of mesh (useful to auto set axes limits)
        ax.scatter(*zip(*self.mesh.points), c="black", s=4, zorder=5)
        num_steps = len(self.states) - 1
        # LEFT MC Plot
        ax_mc.plot(self.phi_list, self.M_list, "k")
        ax.axis("equal")
        # Location indicator as patch so we can update it easily
        # Ellipse b/c we don't have equal axis ranges for phi/M
        point = Ellipse(
            (0, 0), max(self.phi_list) / 50, max(self.M_list) / 50, fc="red", zorder=10
        )
        ax_mc.add_patch(point)
        # Scientific notation for phi values
        ax_mc.ticklabel_format(
            axis="x", style="sci", scilimits=(-2, 2), useMathText=True
        )

        # Setup Sliders for intractability to look through results
        def update(state: float) -> None:
            state_id = int(state)
            self.state_id = int(state_id)
            # Update Patch  colors
            state_colors = [
                self.color_from_state(state)
                for state in self.states[state_id].mat_state
            ]
            for patch, color in zip(patches, state_colors):
                patch.set_facecolor(color)
            point.center = (self.phi_list[state_id], self.M_list[state_id])
            # fig.canvas.draw_idle()

        ax_slider = plt.axes((0.117, 0.01, 0.79, 0.02), facecolor="White")
        slider_fct = Slider(
            ax_slider, "STEP", 0, num_steps, valinit=0, valstep=1, valfmt="%i"
        )
        slider_fct.on_changed(update)

        def frame(state: float) -> list:
            slider_fct.set_val(state)
            return []

        print(
            "Saving animated gif. This will take a while and take several dozen megs of space.\n"
            "It's saving full frames and not optimizing."
        )
        test = anim.FuncAnimation(fig, frame, frames=len(self.states))
        test.save("Moment_Curvature.gif", fps=16)

    def display_mc(self) -> None:
        """Plot the moment curvature, and interactive section viewer"""
        # TODO SEPARATE PLOT FOR FORCE/STRESS DISTRIBUTION?
        # TODO PATCH COLLECTION FOR INSTANTANEOUS? COLOR UPDATE?
        fig, (ax, ax_mc) = plt.subplots(1, 2, figsize=(self.figsize * 2, self.figsize))
        polygons = []
        for element in self.mesh.elements:
            poly_pts = [self.mesh.points[i] for i in element]
            polygons.append(poly_pts)
        patches: list[Patch] = []
        # Draw Polygons for mesh and reinforcement
        i = 0
        for p in polygons:
            patch = plt.Polygon(
                tuple(p), fc="Gray", ec="black", zorder=0, picker=0.01, gid=str(i)
            )
            patches.append(patch)
            ax.add_patch(patch)
            i += 1
        if self.reinforcement:
            for reinforce_group in self.reinforcement:
                for centroid in reinforce_group.points:
                    new_patch = CirclePolygon(centroid, reinforce_group.radius)
                    new_patch.set_facecolor("White")
                    new_patch.set_edgecolor("Black")
                    new_patch.set_picker(0.01)
                    new_patch.set_gid(str(i))
                    new_patch.set_zorder(15)
                    patches.append(new_patch)
                    ax.add_patch(new_patch)
                    i += 1
        # Plot nodal points of mesh (useful to auto set axes limits)
        ax.scatter(*zip(*self.mesh.points), c="black", s=4, zorder=5)
        num_steps = len(self.states) - 1
        # LEFT MC Plot
        ax_mc.plot(self.phi_list, self.M_list, "k")
        ax.axis("equal")
        # Location indicator as patch so we can update it easily
        # Ellipse b/c we don't have equal axis ranges for phi/M
        point = Ellipse(
            (0, 0), max(self.phi_list) / 50, max(self.M_list) / 50, fc="red", zorder=10
        )
        ax_mc.add_patch(point)
        # Scientific notation for phi values
        ax_mc.ticklabel_format(
            axis="x", style="sci", scilimits=(-2, 2), useMathText=True
        )

        # Setup Sliders for intractability to look through results
        def update(phi_step: float) -> None:
            phi_step = int(phi_step)
            self.state_id = int(phi_step)
            # Update Patch  colors
            state_colors = [
                self.color_from_state(state)
                for state in self.states[phi_step].mat_state
            ]
            for patch, color in zip(patches, state_colors):
                patch.set_facecolor(color)
            point.center = (self.phi_list[phi_step], self.M_list[phi_step])
            fig.canvas.draw_idle()

        ax_slider = plt.axes((0.117, 0.01, 0.79, 0.02), facecolor="White")
        slider_fct = Slider(
            ax_slider, "STEP", 0, num_steps, valinit=0, valstep=1, valfmt="%i"
        )
        slider_fct.on_changed(update)

        # Setup mesh onclick event to view stress/strain location of each fiber
        def onclick(event: Event) -> bool:
            if type(event) != PickEvent:
                return False
            gid = event.artist.get_gid()
            patch_id = int(gid if gid else 0)
            strain = self.states[self.state_id].strains[patch_id]
            stress = self.states[self.state_id].stresses[patch_id]
            self.display_material(self.ele_mat[patch_id], (strain, stress))
            return True

        fig.canvas.mpl_connect("pick_event", onclick)

        plt.show()

    def display_mc_2x2(self) -> None:
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
        patches: list[Patch] = []
        # Draw Polygons for mesh and reinforcement
        i = 0
        for p in polygons:
            patch = CirclePolygon(
                tuple(p), fc="Gray", ec="black", zorder=0, picker=0.01, gid=str(i)
            )
            patches.append(patch)
            ax.add_patch(patch)
            i += 1
        if self.reinforcement:
            for reinforce_group in self.reinforcement:
                for centroid in reinforce_group.points:
                    new_patch = CirclePolygon(centroid, reinforce_group.radius)
                    new_patch.set_facecolor("White")
                    new_patch.set_edgecolor("Black")
                    new_patch.set_picker(0.01)
                    new_patch.set_gid(str(i))
                    new_patch.set_zorder(15)
                    patches.append(new_patch)
                    ax.add_patch(new_patch)
                    i += 1
        # Plot nodal points of mesh (useful to auto set axes limits)
        ax.scatter(*zip(*self.mesh.points), c="black", s=4, zorder=5)
        num_steps = len(self.states) - 1
        # LEFT MC Plot
        ax_mc.plot(self.phi_list, self.M_list, "k")
        ax.axis("equal")
        # Location indicator as patch so we can update it easily
        # Ellipse b/c we don't have equal axis ranges for phi/M
        point = Ellipse(
            (0, 0), max(self.phi_list) / 50, max(self.M_list) / 50, fc="red", zorder=10
        )
        ax_mc.add_patch(point)
        # Scientific notation for phi values
        ax_mc.ticklabel_format(
            axis="x", style="sci", scilimits=(-2, 2), useMathText=True
        )

        line = line2 = line3 = None

        # Setup Sliders for intractability to look through results
        def update(phi_step: float) -> None:
            phi_step = int(phi_step)
            self.state_id = int(phi_step)
            # Update Patch  colors
            state_colors = [state.name for state in self.states[phi_step].mat_state]
            for patch, color in zip(patches, state_colors):
                patch.set_facecolor(color)
            point.center = (self.phi_list[phi_step], self.M_list[phi_step])

            x = [self.states[phi_step].min_strain, self.states[phi_step].max_strain]
            y = [self.analysis_model.min_y, self.analysis_model.max_y]

            if line:
                line.set_data(x, y)
            if line2:
                line2.set_data(
                    [-1, 1], [self.states[phi_step].y_loc, self.states[phi_step].y_loc]
                )
            if line3:
                line3.set_data(
                    [self.analysis_model.min_y * 1.1, self.analysis_model.max_y * 1.1],
                    [self.states[phi_step].y_loc, self.states[phi_step].y_loc],
                )
            fig.canvas.draw_idle()

        # Setup mesh onclick event to view stress/strain location of each fiber
        def onclick(event: Event) -> bool:
            if type(event) != PickEvent:
                return False
            gid = event.artist.get_gid()
            patch_id = int(gid if gid else 0)
            strain = self.states[self.state_id].strains[patch_id]
            stress = self.states[self.state_id].stresses[patch_id]
            self.display_material(self.ele_mat[patch_id], (strain, stress))
            return True

        fig.canvas.mpl_connect("pick_event", onclick)

        ax_slider = plt.axes((0.117, 0.01, 0.79, 0.02), facecolor="white")
        slider_fct = Slider(
            ax_slider, "STEP", 0, num_steps, valinit=0, valstep=1, valfmt="%i"
        )
        slider_fct.on_changed(update)

        # Strain plot
        ax_strain.plot(
            [0, 0], [self.analysis_model.min_y, self.analysis_model.max_y], "k"
        )
        ax_strain.set_xlim(
            [
                self.states[len(self.states) - 1].min_strain * 1.05,
                self.states[len(self.states) - 1].max_strain * 1.05,
            ]
        )

        (line,) = ax_strain.plot(
            [0, 0], [self.analysis_model.min_y, self.analysis_model.max_y], color="Red"
        )
        phi_step_init = 0
        (line2,) = ax_strain.plot(
            [-1, 1],
            [self.states[phi_step_init].y_loc, self.states[phi_step_init].y_loc],
            color="Blue",
        )
        (line3,) = ax.plot(
            [self.analysis_model.min_y * 1.1, self.analysis_model.max_y * 1.1],
            [self.states[phi_step_init].y_loc, self.states[phi_step_init].y_loc],
            color="Red",
        )

        # Stress Plot
        ax_stress.plot(
            [0, 0], [self.analysis_model.min_y, self.analysis_model.max_y], "k"
        )
        ax_stress.text(0, 0, "This Plot Was\nIntentionally Left Blank\n(For Now)")

        plt.tight_layout()
        plt.show()

    def export_results(self, filename: str = "") -> None:
        out = ["phi,M"]
        for phi, M in zip(self.phi_list, self.M_list):
            out.append(f"{phi},{M}")
        if filename == "":
            print("\n".join(out))
        else:
            with open(filename, "w") as f:
                f.write("\n".join(out))
