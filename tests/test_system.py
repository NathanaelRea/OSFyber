import sys
import unittest
from io import StringIO

from osfyber.materials import (
    ConfinedConcreteProps,
    ConfinedPressureCircleProps,
    conf_pressure_circle,
)
from osfyber.system import (
    CircleGeometryProps,
    FyberModel,
    LoadProps,
    ReinforcementProps,
    SteelProps,
    UnconfinedConcreteProps,
)


class Capturing(list):
    def __enter__(self) -> "Capturing":
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, exc_type: type, exc_value: Exception, traceback: object) -> None:
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio
        sys.stdout = self._stdout


class TestFyberModel(unittest.TestCase):
    def test_analyze(self) -> None:
        model = FyberModel()
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
        model.add_geometry_circle(
            unconfined_concrete, CircleGeometryProps(D=36, c=(0, 0))
        )
        model.add_reinforcement_circle(
            steel,
            CircleGeometryProps(D=29.5, c=(0, 0)),
            ReinforcementProps(bar=9, count=12, conf_material=confined_concrete),
        )
        model.set_load(LoadProps(P=1000))
        model.generate_mesh()

        with Capturing() as output:
            model.analyze()

        expected = [
            "Analysis ended at phi=0.00166",
            "Failure: Confined Concrete Crushing",
            "\tMax Available Strain=-0.01432",
            "\tStrain Experienced=-0.01443",
            "\tMat_id=2",
            "\tLocation=[-1.148, 13.95]",
        ]

        self.assertEqual(expected, output)


if __name__ == "__main__":
    unittest.main()
