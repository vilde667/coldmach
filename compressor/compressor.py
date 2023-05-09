class Compressor:

    mark: str
    diameter: float
    piston_stroke: float
    theor_v: float
    frequency: int
    Q_R717: float
    N_R717: float

    def __init__(self, params: dict) -> None:
        self.mark = params['mark']
        self.diameter = params['diameter']
        self.piston_stroke = params['piston_stroke']
        self.theor_v = params['theor_v']
        self.frequency = params['frequency']
        self.Q_R717 = params['Q_R717']
        self.N_R717 = params['N_R717']
