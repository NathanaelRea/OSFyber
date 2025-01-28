# OSFyber

Open Source 'Fyber' Analysis Tool for Moment Curvature.

![Example 1 MC Scrub](Pics/Example_1_Scrub.gif)

## Installation

- uv init
- uv add git+<https://github.com/NathanaelRea/OSFyber>

## Examples

### Confined Circular Column

[circular_column.py](/tests/circular_column.py)

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

You can click on the mesh to view stress strain plot at that step:

![Example click on mesh](Pics/Example_click.png)
