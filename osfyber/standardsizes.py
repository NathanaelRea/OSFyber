"""
A place to store standard structural properties, such as:
    Reinforcement Sizes and Areas
    Concrete and Steel yield strength
"""


class Standard:
    def bar(self, bar_number):
        bars = {3: (0.375, 0.110), 4: (0.500, 0.200), 5: (0.625, 0.310),
                6: (0.750, 0.440), 7: (0.875, 0.600), 8: (1.000, 0.790), 9: (1.128, 1.000),
                10: (1.270, 1.270), 11: (1.410, 1.560), 14: (1.693, 2.250), 18: (2.257, 4.000)}
        return Bar(bars[bar_number][0], bars[bar_number][1])


class Bar:
    def __init__(self, db, area):
        self.db = db
        self.area = area
