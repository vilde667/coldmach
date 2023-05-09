import numpy as np
import math
import sqlite3

from connector import DbConnector
c = DbConnector('./db/main.db').c

import numpy as np
import pandas as pd
from scipy import signal
from scipy import interpolate

from params import Param
from params import ParamFormula

from compressor.compressor import Compressor
from condenser.condenser import Condenser


class MainSolver:
    def __init__(self, Q0, city, agent, temp_brine_out, condenser_type, water_supply_type) -> None:
        self.Q0 = Param('Q0', 'Расчетная холодопроизводительность')
        self.Q0.set_value(Q0, 'кВт')

        self.city = city

        self.agent = agent

        self.t_brine_out = Param('t\'\'р', 'Температура рассола на выходе')
        self.t_brine_out.set_value(temp_brine_out, '°С')

        self.condenser_type = condenser_type

        self.water_supply_type = water_supply_type

        self.set_delta_t(5.0)
        self.set_delta_t_p(3.0)
        self.set_t_brine_in()
        self.set_delta_t_v(3.0)
        self.set_t_in_v(self.city)
        self.set_t_out_v()
        self.set_t_0()
        self.set_delta_t_k(5.0)
        self.set_t_k()
        self.set_delta_t_p(4.0)
        self.set_t_p()
        self.set_delta_t_vs(8.0)
        self.set_t_vs()

        self.graphic = Graphic(self.t_0.value, self.t_vs.value, self.t_k.value, self.t_p.value)

        self.graphic.make_graphic_template()
        self.graphic.make_graphic()

        self.set_q_0()
        self.set_l()
        self.set_q_k()
        self.set_e()
        self.set_exerg()
        self.set_G_0()
        self.set_V_0()
        self.set_q_v()

        self.standart_graphic = Graphic(-15, -10, 30, 25)

        self.standart_graphic.make_graphic_template()
        self.standart_graphic.make_graphic()

        self.set_q_0_st()
        self.set_q_v_st()
        self.set_compression_ratio()
        self.set_compression_ratio_standart()

        self.set_lambda_real()
        self.set_lambda_standart()
        self.set_Q0_st()

        self.set_compressor()

        self.set_N()
        self.set_n()
        self.set_G_0_()

        self.condenser = Condenser(self.compressor, self.agent, self.delta_t_v, self.t_in_v, self.t_out_v, self.t_k)

    def set_delta_t(self, value):
        self.delta_t = Param('Δt0', 'разность температур (средняя температура рассола - кипение хладагента)')
        self.delta_t.set_value(value, '°С')
        self.delta_t.set_comment('Принимаем по таблице 4.1 пособия')
        self.delta_t.use()

    def set_delta_t_p(self, value):
        self.delta_t_p = Param('Δtр', 'разность температур (рассол на входе - рассол на выходе)')
        self.delta_t_p.set_value(value, '°С')
        self.delta_t_p.set_comment('Принимаем по таблице 4.1 пособия')
        self.delta_t_p.use()

    def set_t_brine_in(self):
        self.t_brine_in = ParamFormula('t\'р', 'температура рассола на входе в ')
        self.t_brine_in.set_formula('delta_t_p + t_brine_out', {'delta_t_p': self.delta_t_p.value, 't_brine_out': self.t_brine_out.value})
        self.t_brine_in.set_unit('°С')
        self.t_brine_in.use()

    def set_delta_t_v(self, value):
        self.delta_t_v = Param('Δtв', 'разность температур (выходящая вода конденсатора - входящая вода конденсатора)')
        self.delta_t_v.set_value(value, '°С')
        self.delta_t_v.set_comment('Принимаем по таблице 4.1 пособия')
        self.delta_t_v.use()

    def set_t_in_v(self, city):

        c.execute('SELECT t_air_out FROM citys WHERE city=\'{0}\''.format(city))
        res = c.fetchall()
        if len(res) == 1:
            result = res[0][0]
        else:
            result = None


        self.t_in_v = Param('t\'в', 'температура воды на входе в конденсатор')
        self.t_in_v.set_value(result, '°С')
        self.t_in_v.set_comment('Находим по темепратуре наружного воздуха в теплый период года по городу установки оборудования - {0}'.format(city))
        self.t_in_v.use()

    def set_t_out_v(self):
        self.t_out_v = ParamFormula('t\'\'в', 'температура воды на выходе из конденсатора')
        self.t_out_v.set_formula('t_in_v + delta_t_v', {'t_in_v': self.t_in_v.value, 'delta_t_v': self.delta_t_v.value})
        self.t_out_v.set_unit('°С')
        self.t_out_v.use()

    def set_t_0(self):
        self.t_0 = ParamFormula('t0', 'температура кипения хладагента в испарителе')
        self.t_0.set_formula('(t_brine_in + t_brine_out) / 2 - delta_t', {'t_brine_in': self.t_brine_in.value, 't_brine_out': self.t_brine_out.value, 'delta_t': self.delta_t.value})
        self.t_0.set_unit('°С')
        self.t_0.use()

    def set_delta_t_k(self, value):
        self.delta_t_k = Param('Δtk', 'разность температур (конденсация хладагента - охлаждающая вода)')
        self.delta_t_k.set_value(value, '°С')
        self.delta_t_k.set_comment('Принимаем по таблице 4.1 пособия')
        self.delta_t_k.use()

    def set_t_k(self):
        self.t_k = ParamFormula('tk', 'температура конденсации')
        self.t_k.set_formula('(t_out_v + t_in_v) / 2 + delta_t_k', {'delta_t_k': self.delta_t_k.value, 't_out_v': self.t_out_v.value, 't_in_v': self.t_in_v.value})
        self.t_k.set_unit('°С')
        self.t_k.use()

    def set_delta_t_p(self, value):
        self.delta_t_p = Param('Δtk', 'разность температур (темепература конденсации - температура переохлаждения)')
        self.delta_t_p.set_value(value, '°С')
        self.delta_t_p.set_comment('Принимаем по таблице 4.1 пособия')
        self.delta_t_p.use()

    def set_t_p(self):
        self.t_p = ParamFormula('tp', 'температура переохлаждения')
        self.t_p.set_formula('t_k - delta_t_p', {'t_k': self.t_k.value, 'delta_t_p': self.delta_t_p.value})
        self.t_p.set_unit('°С')
        self.t_p.use()

    def set_delta_t_vs(self, value):
        self.delta_t_vs = Param('Δtвс', 'разность температур (хладагент засасываемый компрессором - кипение хладагента)')
        self.delta_t_vs.set_value(value, '°С')
        self.delta_t_vs.set_comment('Принимаем по таблице 4.1 пособия')
        self.delta_t_vs.use()

    def set_t_vs(self):
        self.t_vs = ParamFormula('tвс', 'температура всасывания')
        self.t_vs.set_formula('t_0 + delta_t_vs', {'t_0': self.t_0.value, 'delta_t_vs': self.delta_t_vs.value})
        self.t_vs.set_unit('°С')
        self.t_vs.use()

    def set_q_0(self):
        self.q_0 = ParamFormula('q0', 'удельная холодопроизводительность')
        self.q_0.set_formula('h1_ - h4', {'h1_': self.graphic.dot_1_.H, 'h4': self.graphic.dot_4.H})
        self.q_0.set_unit('кДж/кг')
        self.q_0.use()

    def set_l(self):
        self.l = ParamFormula('l', 'удельная работа сжатия в компрессоре')
        self.l.set_formula('h2 - h1', {'h2': self.graphic.dot_2.H, 'h1': self.graphic.dot_1.H})
        self.l.set_unit('кДж/кг')
        self.l.use()

    def set_q_k(self):
        self.q_k = ParamFormula('qk', 'удельная тепловая нагрузка на конденсатор')
        self.q_k.set_formula('h2 - h3_', {'h2': self.graphic.dot_2.H, 'h3_': self.graphic.dot_3_.H})
        self.q_k.set_unit('кДж/кг')
        self.q_k.use()

    def set_e(self):
        self.e = ParamFormula('e', 'холодильный коэффициент цикла')
        self.e.set_formula('(h1_-h4)/(h2-h1)', {'h1_': self.graphic.dot_1_.H, 'h4': self.graphic.dot_4.H, 'h2': self.graphic.dot_2.H, 'h1': self.graphic.dot_1.H})
        self.e.set_unit('')
        self.e.use()

    def set_exerg(self):
        self.exerg = ParamFormula('ηe', 'эксергетический КПД')
        self.exerg.set_formula('(T0 - Tos)/T0 * (q0/l)', {'T0': self.graphic.dot_4.T + 273.15, 'Tos': self.t_in_v.value + 273.15, 'q0': self.q_0.value, 'l': self.l.value})
        self.exerg.set_unit('')
        self.exerg.use()

    def set_G_0(self):
        self.G_0 = ParamFormula('G0', 'масса циркулирующего хладагента')
        self.G_0.set_formula('Q0/q0', {'Q0': self.Q0.value, 'q0': self.q_0.value})
        self.G_0.set_unit('кг/с')
        self.G_0.use()

    def set_V_0(self):
        self.V_0 = ParamFormula('V0', 'Действительный объем пара')
        self.V_0.set_formula('G0*teta', {'G0': self.G_0.value, 'teta': self.graphic.dot_1.V})
        self.V_0.set_unit('м3/с')
        self.V_0.use()

    def set_q_v(self):
        self.q_v = ParamFormula('qv', 'объемная холодопроизводительность')
        self.q_v.set_formula('q0/teta', {'q0': self.q_0.value, 'teta': self.graphic.dot_1.V})
        self.q_v.set_unit('кДж/м3')
        self.q_v.use()

    def set_q_0_st(self):
        self.q_0_st = ParamFormula('q0 ст', 'удельная холодопроизводительность стандартного режима')
        self.q_0_st.set_formula('h1_ - h4', {'h1_': self.standart_graphic.dot_1_.H, 'h4': self.standart_graphic.dot_4.H})
        self.q_0_st.set_unit('кДж/кг')
        self.q_0_st.use()

    def set_q_v_st(self):
        self.q_v_st = ParamFormula('qv ст', 'объемная холодопроизводительность стандартного режима')
        self.q_v_st.set_formula('q0/teta', {'q0': self.q_0_st.value, 'teta': self.standart_graphic.dot_1.V})
        self.q_v_st.set_unit('кДж/м3')
        self.q_v_st.use()

    def set_compression_ratio(self):
        self.compression_ratio = ParamFormula('Pk/P0', 'степень сжатия действительного компрессора')
        self.compression_ratio.set_formula('Pk/P0', {'Pk': self.graphic.dot_2.P, 'P0': self.graphic.dot_1.P})
        self.compression_ratio.set_unit('')
        self.compression_ratio.use()

    def set_compression_ratio_standart(self):
        self.compression_ratio_standart = ParamFormula('Pk/P0 ст', 'степень сжатия стандартного компрессора')
        self.compression_ratio_standart.set_formula('Pk/P0', {'Pk': self.standart_graphic.dot_2.P, 'P0': self.standart_graphic.dot_1.P})
        self.compression_ratio_standart.set_unit('')
        self.compression_ratio_standart.use()

    def set_lambda_real(self):
        value = TableValue(c, 'R717', {'compression_ratio': self.compression_ratio.value}, 'feed_ratio')
        self.lambda_real = Param('λ', 'коэффициент подачи реального компрессора')
        self.lambda_real.set_value(value.get_value(), '')
        self.lambda_real.set_comment('Находим по графику 4.11 для R717 по степени сжатия={0}'.format(self.compression_ratio.value))
        self.lambda_real.use()

    def set_lambda_standart(self):
        value = TableValue(c, 'R717', {'compression_ratio': self.compression_ratio_standart.value}, 'feed_ratio')
        self.lambda_standart = Param('λ', 'коэффициент подачи стандартного компрессора')
        self.lambda_standart.set_value(value.get_value(), '')
        self.lambda_standart.set_comment('Находим по графику 4.11 для R717 по степени сжатия={0}'.format(self.compression_ratio_standart.value))
        self.lambda_standart.use()

    def set_Q0_st(self):
        self.Q0_st = ParamFormula('Q0 ст', 'стандартная холодопроизводительность компрессора')
        self.Q0_st.set_formula('Q0 * (q_v_st * lambda_standart)/(q_v * lambda_real)', {'Q0': self.Q0.value, 'q_v_st': self.q_v_st.value, 'q_v': self.q_v.value, 'lambda_standart': self.lambda_standart.value, 'lambda_real': self.lambda_real.value})
        self.Q0_st.set_unit('кВт')
        self.Q0_st.use()

    def set_compressor(self):
        c.execute('SELECT mark,diameter,piston_stroke,theor_v,frequency,Q_R717,N_R717 FROM compressors WHERE Q_R717 IS NOT NULL AND Q_R717>{0} ORDER BY Q_R717'.format(self.Q0_st.value))
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

            self.compressor = Compressor(params)

            print('По полученной холодопроизводительности выбираем компрессор марки ' + self.compressor.mark)
    
    def set_N(self):
        self.N = Param('N', 'количество компрессоров холодильной установки')
        self.N.set_value(math.ceil(self.Q0_st.value/self.compressor.Q_R717), '')
        self.N.set_comment('Расчитываем')
        self.N.use()

    def set_n(self):
        self.n = Param('n', 'количество компрессоров холодильной установки с учетом резерва')
        self.n.set_value(self.N.value + 1, '')
        self.n.set_comment('Расчитываем')
        self.n.use()

    def set_G_0_(self):
        self.G_0_ = ParamFormula('G\'0', 'масса циркулирующего хладагента на один компрессор')
        self.G_0_.set_formula('G0/n', {'G0': self.G_0.value, 'n': self.n.value})
        self.G_0_.set_unit('кг/с')
        self.G_0_.use()





class TableValue:
    def __init__(self, cursor:sqlite3.Cursor, param:dict, known_param:dict, table_name:str) -> None:
        self.param = param
        self.known_param_key = list(known_param.keys())[0]
        self.known_param_value = known_param[self.known_param_key]
        self.table_name = table_name
        self.cursor = cursor

    def get_value(self):
        known_values = []
        c.execute('SELECT {0} FROM {1}'.format(self.known_param_key, self.table_name))
        res = c.fetchall()
        for r in res:
            known_values.append(float(r[0]))
        res = list(set(known_values))
        res.sort()

        known_value_min = min(res, key=lambda x:abs(x-self.known_param_value))

        if known_value_min > self.known_param_value:
            known_value_min = res[res.index(known_value_min) - 1]

        known_value_max = res[res.index(known_value_min) + 1]

        if self.known_param_value == known_value_min:
            c.execute('SELECT {0} FROM {1} WHERE {2}={3}'.format(self.param, self.table_name, self.known_param_key, known_value_min))
            res = c.fetchall()
            for r in res:
                return r[0]
        elif self.known_param_value == known_value_max:
            c.execute('SELECT {0} FROM {1} WHERE {2}={3}'.format(self.param, self.table_name, self.known_param_key, known_value_max))
            res = c.fetchall()
            for r in res:
                return r[0]
        else:
            param_min = 0
            param_max = 0
            c.execute('SELECT {0} FROM {1} WHERE {2}={3}'.format(self.param, self.table_name, self.known_param_key, known_value_min))
            res = c.fetchall()
            for r in res:
                param_min = r[0]

            c.execute('SELECT {0} FROM {1} WHERE {2}={3}'.format(self.param, self.table_name, self.known_param_key, known_value_max))
            res = c.fetchall()
            for r in res:
                param_max = r[0]

            k = (self.known_param_value - known_value_min) / (known_value_max - known_value_min)
            return param_min + (param_max - param_min)*k



        


        

solver = MainSolver(180.0, 'Владимир', 'R717', -2.0, 'КТ', 'оборотное')

    
