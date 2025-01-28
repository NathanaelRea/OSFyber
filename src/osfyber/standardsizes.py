class Bar:
    def __init__(self, db: float, area: float):
        self.db = db
        self.area = area


BARS = {
    3: Bar(0.375, 0.110),
    4: Bar(0.500, 0.200),
    5: Bar(0.625, 0.310),
    6: Bar(0.750, 0.440),
    7: Bar(0.875, 0.600),
    8: Bar(1.000, 0.790),
    9: Bar(1.128, 1.000),
    10: Bar(1.270, 1.270),
    11: Bar(1.410, 1.560),
    14: Bar(1.693, 2.250),
    18: Bar(2.257, 4.000),
}


class Standard:
    def bar(self, bar_number: int) -> Bar:
        if bar_number not in BARS:
            raise ValueError(f"Bar number {bar_number} does not exist")
        return BARS[bar_number]
