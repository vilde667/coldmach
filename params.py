import numexpr as ne


class Param:
    name: str
    desc: str
    value: float
    unit: str
    comment: str

    def __init__(self, name: str, desc: str) -> None:
        self.name = name
        self.desc = desc

    def set_value(self, value: float, unit: str):
        self.value = value
        self.unit = unit

    def set_comment(self, comment: str):
        self.comment = comment

    def use(self):
        print('{0} = {1} {2} ({3})'.format(self.name, self.value, self.unit, self.comment))


class ParamFormula:
    name: str
    desc: str
    value: float
    unit: str
    comment: str
    formula: str
    params: dict

    def __init__(self, name: str, desc: str) -> None:
        self.name = name
        self.desc = desc

    def set_unit(self, unit: str):
        self.unit = unit

    def set_comment(self, comment: str):
        self.comment = comment

    def set_formula(self, formula: str, params: dict):
        self.formula = formula
        self.params = params

    def use(self):
        self.value = ne.evaluate(self.formula, local_dict=self.params)
        self.comment = 'Вычисляем по формуле {0}'.format(self.formula)
        print(self.desc)
        print('{0} = {1} {2} ({3})'.format(self.name, self.value, self.unit, self.comment))
