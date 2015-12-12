"""
Classes which parse marshalled data. Some parsers read data from files. Others
read over the network. All parsers share the same interface.
"""

class Parser(object):
    def __init__(self):
        self.n_ = 0       # number of time steps in particle simulation
        self.L_ = 0.0     # side length of particle space
        self.N_ = 0       # side length of particle grid/mesh
        self.N_p_ = 0     # number of particles
        self.masses_ = [] # masses of particles
        self.last_ = []   # the last set of data read

    @property
    def n(self):
        return self.n_

    @property
    def L(self):
        return self.L_

    @property
    def N(self):
        return self.N_

    @property
    def N_p(self):
        return self.N_p_

    @property
    def masses(self):
        return self.masses_

    def next(self):
        """
        Parses the next set of positions. If the end of the stream is reached,
        then the last valid data is returned.
        """
        raise NotImplementedError("next not implemented in interface")

class FileParser(Parser):
    def __init__(self, filename):
        self.f_ = open(filename, "r")
        self.n_ = int(self.f_.readline())
        self.L_ = float(self.f_.readline())
        self.N_ = int(self.f_.readline())
        self.N_p_ = int(self.f_.readline())
        self.masses_ = [float(self.f_.readline()) for _ in range(self.N_p)]
        self.last_ = None

    def next(self):
        try:
            positions = []
            for _ in range(self.N_p):
                x, y, _, _ = [float(_) for _ in self.f_.readline().split()]
                positions.append((x, y))
            self.last_ = positions
            return self.last_
        except ValueError:
            # A ValueError is raised when we try to read off the end of the
            # file. We catch this exception and return the last valid data
            # read.
            return self.last_


class NetParser(Parser):
    def __init__(self, host, port):
        raise NotImplementedError("NetParser not yet implemented")

    def next(self):
        raise NotImplementedError("NetParser not yet implemented")
