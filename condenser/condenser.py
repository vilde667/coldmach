from params import Param
from params import ParamFormula

from compressor.compressor import Compressor


class Condenser:
    
    compressor: Compressor
    coolant: str
    K: Param
    Qk: ParamFormula
    Cp: Param
    tv_in: Param
    tv_out: ParamFormula
    tk: ParamFormula
    tv_delta: Param
    delta_t_log: ParamFormula

    def __init__(self, compressor: Compressor, coolant: str, delta_t_v: Param, t_in_v: Param, t_out_v: ParamFormula,
                 t_k: ParamFormula) -> None:
        self.compressor = compressor
        self.coolant = coolant
        self.tv_delta = delta_t_v
        self.tv_in = t_in_v
        self.tv_out = t_out_v
        self.tk = t_k

        self.set_K()
        self.set_Qk()
        self.set_Cp()
        self.set_Gv()
        self.set_delta_t_min_max()
        self.set_delta_t_log()
        self.set_F_prev()

    def set_K(self):
        self.K = Param('K', 'коэффициент теплопередачи заданный')
        if self.coolant == 'R717':
            self.K.set_value(750, 'Вт/(м2*°С)')
        else:
            self.K.set_value(2100, 'Вт/(м2*°С)')
        self.K.set_comment('Выбираем по таблице 4.14')
        self.K.use()

    def set_Qk(self):
        self.Qk = ParamFormula('Qk', 'Нагрузка на конденсатор')
        self.Qk.set_formula('Q0comp + Ncomp', {'Q0comp': self.compressor.Q_R717, 'Ncomp': self.compressor.N_R717})
        self.Qk.set_unit('кВт')
        self.Qk.use()

    def set_Cp(self):
        # Сделать расчетным по таблице!
        self.Cp = Param('Cp', 'Удельная теплоемкость воды в конденсаторе')
        self.Cp.set_value(4200, 'Дж/(кг*°С)')
        self.Cp.set_comment('Приблизительное, при температуре 5°С')
        self.Cp.use()

    def set_Gv(self):
        self.Gv = ParamFormula('Gv', 'Расход воды на конденсатор')
        self.Gv.set_formula('Qk/(Cp*delta_tv)', {'Qk': self.Qk.value*1000, 'Cp': self.Cp.value, 'delta_tv': self.tv_delta.value})
        self.Gv.set_unit('кг/с')
        self.Gv.use()

    def set_delta_t_min_max(self):
        deltas = []
        deltas.append(self.tk.value - self.tv_in.value)
        deltas.append(self.tk.value - self.tv_out.value)

        self.tmin_delta = min(deltas)
        self.tmax_delta = max(deltas)

    def set_delta_t_log(self):
        self.delta_t_log = ParamFormula('Δtlog', 'среднелогарифмическая разность температур по воде')
        self.delta_t_log.set_formula('(tmax_delta - tmin_delta)/(log(tmax_delta/tmin_delta))', {'tmax_delta': self.tmax_delta, 'tmin_delta': self.tmin_delta})
        self.delta_t_log.set_unit('°С')
        self.delta_t_log.use()

    def set_F_prev(self):
        self.F_prev = ParamFormula('F', 'предварительно расчитаная площадь конденсатора')
        self.F_prev.set_formula('Qk/(K * delta_t_log)', {'Qk': self.Qk.value*1000, 'K': self.K.value, 'delta_t_log': self.delta_t_log.value})
        self.F_prev.set_unit('м2')
        self.F_prev.use()

    def set_condensator_by_table(self):
        c.execute('SELECT mark,S,diam_out,lenght,width,height,pipes,pipes_lenght,strokes,pipes_diam_in,pipes_diam_out FROM condenser_horisontal_amm WHERE S IS NOT NULL AND S>{0} ORDER BY Q_R717'.format(self.F_prev.value))
        row = c.fetchone()
        if row:
            params = {}
            params['mark'] = row[0]
            params['diameter'] = row[1]
            params['piston_stroke'] = row[2]
            params['theor_v'] = row[3]
            params['frequency'] = row[4]
            params['Q_R717'] = row[5]
            params['N_R717'] = row[6]

            print('По полученной площади поверхности выбираем конденсатор марки ' + self.compressor.mark)

    def set_A(self):
        self.delta_t_log = ParamFormula('A', 'Числовой коэффициент для расчета удельного теплового потока')
        self.delta_t_log.set_formula('0.728 * B * d_out**-0.25 * n_avg**-0.167', {'B': self.tmax_delta, 'd_out': self.tmin_delta, 'n_avg': self.tmin_delta})
        self.delta_t_log.set_unit('°С')
        self.delta_t_log.use()