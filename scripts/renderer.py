import numpy as np

def render(positions, masses, max_mass, a):
    # TODO: handle masses
    for (x, y) in positions:
        try:
            a[x, y] = 255
        except IndexError as e:
            # TODO: make this impossible
            print e

class Renderer(object):
    def __init__(self, pixels, max_mass, parser):
        self.pixels_ = pixels
        self.parser_ = parser
        self.max_mass_ = max_mass
        self.a_ = np.zeros((self.pixels_, self.pixels_, 3))

    @property
    def pixels(self):
        return self.pixels_

    def scale(self, x):
        return int(x * float(self.pixels_) / self.parser_.L)

    def next(self):
        positions = self.parser_.next()
        positions = [(self.scale(x), self.scale(y)) for (x, y) in positions]
        self.a_[:] = 0
        render(positions, self.parser_.masses, self.max_mass_, self.a_)
        return self.a_
