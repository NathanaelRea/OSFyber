from osfyber.materials import (
    ConfinedConcreteProps,
    ConfinedPressureCircleProps,
    SteelProps,
    UnconfinedConcreteProps,
    conf_pressure_circle,
)
from osfyber.system import (
    CircleGeometryProps,
    FyberModel,
    LoadProps,
    ReinforcementProps,
)

model = FyberModel()

# Set Materials
unconfined_concrete = model.add_material_concrete(
    UnconfinedConcreteProps(fpc=5.2, ecp=-0.002, ecu=-0.005)
)
fple = conf_pressure_circle(
    ConfinedPressureCircleProps(fyh=68, bar_number=4, D=31.5, s=3)
)
confined_concrete = model.add_material_confined_concrete(
    ConfinedConcreteProps(fpc=5.2, fple=fple)
)
steel = model.add_material_steel(
    SteelProps(E=29565, fy=68, fsu=95, e_sh=0.0125, e_su=0.09, P=2.8)
)

# Set geometry
model.add_geometry_circle(unconfined_concrete, CircleGeometryProps(D=36, c=(0, 0)))
model.add_reinforcement_circle(
    steel,
    CircleGeometryProps(D=29.5, c=(0, 0)),
    ReinforcementProps(bar=9, count=12, conf_material=confined_concrete),
)

# Add loading (positive is compression)
model.set_load(LoadProps(P=1000))

# Generate mesh
model.generate_mesh()

# Check input
# model.display_materials()
# model.display_mesh()

model.analyze()
model.display_mc()
# model.display_mc_2x2()

# save results
# model.export_results()
