class Point:
    P: float = 0
    T: float = 0
    H: float = 0
    S: float = 0
    V: float = 0

    def set_p(self, p: float):
        self.P = p

    def set_t(self, t: float):
        self.T = t

    def set_h(self, h: float):
        self.H = h

    def set_s(self, s: float):
        self.S = s

    def set_v(self, v: float):
        self.V = v

    def get_params(self):
        return self.P, self.T, self.H, self.S, self.V
    
    def print_params(self):
        print(self.P, self.T, self.H, self.S, self.V)
