# OSFyber

Open Source 'Fyber' Analysis Tool for Moment Curvature.

![Example 1 MC Scrub](Pics/Example_1_Scrub.gif)

## Installation

- Ensure you have MS build tools for MeshPy
  - <https://visualstudio.microsoft.com/visual-cpp-build-tools>
- Initialize Poetry project
  - `poetry init`
- Add directly from repo url
  `poetry add git+https://github.com/NathanaelRea/OSFyber`

## Examples

### Standard Confined Circular Column

```python
from osfyber.system import FyberModel

model = FyberModel()

model.add_material(mat_num=1, mat='concrete', fpc=5.2, ecp=-.002, ecu=-.005)
fple = FyberModel.conf_pressure('circle', fyh=68, bar=4, D=31.5, s=3)
model.add_material(mat_num=2, mat='concrete', fpc=5.2, fple=fple)
model.add_material(mat_num=3, mat='steel', E=29565, fy=68, fsu=95, e_sh=0.0125, e_su=0.09, P=2.8)

model.add_geometry('circle', mat_id=1, c=(0, 0), D=36)

model.add_reinforcement('circle', mat_id=3, D=29.5, c=(0, 0), bar=9, count=12, conf_id=2)

model.set_load('Axial', P=1000)

model.generate_mesh()

# Check input materials
# model.display_materials()
# Check mesh
# model.display_mesh()

# Run analysis and display results
model.analyze()

# Save results
# model.export_results()

# Display Moment Curvature
model.display_mc()
```

Output from display_materials():

![Example 1 Material Unconfined Concrete](Pics/Example_1_Mat_1.png)

![Example 1 Material Confined Concrete](Pics/Example_1_Mat_2.png)

![Example 1 Material Steel](Pics/Example_1_Mat_3.png)

Output from analyze():

```text
Analysis ended at phi=0.0017
Failure: Confined Concrete Crushing
    Max Available Strain=-0.01432
    Strain Experienced=-0.01512
    Mat_id=2
    Location=[3.281, 13.664]
```

Fully Interactive Moment Curvature Diagram. Ability to step through process of degradation. Can also click on a mesh to see internal properties. Can use Mouse on either plot to view (x,y) values in bottom right.

Output from display_mc():

![Example 1 Disp Mc](Pics/Example_1_Disp_MC.png)
