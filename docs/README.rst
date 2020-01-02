# OSFyber

Open Source 'Fyber' Analysis Tool for Moment Curvature.

Currently in pre-alpha. Still implementing features in todo, generalizing code, and fixing bugs.

## Installation

For the actively developed version:

```
$ pip install git+https://github.com/NathanaelRea/OSFyber
```

Or just download the zip from github and use pip install on the extracted folder.

## Examples

### Standard Circular Section

```python
from osfyber.system import FyberModel

### CREATE A MODEL
model = FyberModel()

### Set Materials
model.add_material(id=1, type='concrete', fpc=5.2, ecp=-.002, ecu=-.005)
fple = model.conf_pressure('circle', fyh=68, bar=4, D=31.5, s=3)
model.add_material(id=2, type='concrete', fpc=5.2, fple=fple)
model.add_material(id=3, type='steel', E=29565, fy=68, fsu=95, e_sh=0.0125, e_su=0.09, P=2.8)

### Set Geometry
model.add_geometry('circle', mat_id=1, c=(0,0), D=36)

### Set longitudinal reinforcmenet, and optinally confinement
model.add_reinforcement('circle', mat_id=3, D=29.5, c=(0,0), bar=9, count=12, conf_id=2)

## ADD LOADING (Positive is Compression)
model.add_load('Axial', P=1000)

### GENERATE MESH
model.generate_mesh()

### CHECK INPUT
model.display_mat()
#model.display_mesh()

### ANALYZE MODEL
model.analyze()

#model.export_results()

### DISPLAY MOMENT CURVATURE
model.display_mc()
```

Output from analysis:

```
Analysis ended at phi=0.0017
Confined Concrete Crushing
Max Available Strain=-0.01432
Strain Experienced=-0.01512
Mat_id=2
Location=[3.281, 13.664]
```


![Example 1 Material Unconfined Concrete](Pics/Example_1_Mat_1.png)

![Example 1 Material Confined Concrete](Pics/Example_1_Mat_2.png)

![Example 1 Material Steel](Pics/Example_1_Mat_3.png)

![Example 1 Initial Plot](Pics/Example_1_Pic_1.png)

![Example 1 First point within yield](Pics/Example_1_Pic_2.png)

![Example 1 First point with cover crushing on compression end](Pics/Example_1_Pic_3.png)

![Example 1 Ultimate ductility - Failure of confined concrete](Pics/Example_1_Pic_4.png)


## TODO

[Trello](https://trello.com/b/FFrJVfhk/osfyber)

