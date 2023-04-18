import numpy as np
import sqlite3
import numexpr as ne
import math

import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots

import numpy as np
import pandas as pd
from scipy import signal
from scipy import interpolate

def frange(start, stop, step):
  i = start
  while i < stop:
    yield i
    i += step

conn = sqlite3.connect('main.db')
c = conn.cursor()

class Param:
    name:str
    desc:str
    value:float
    unit:str
    comment:str

    def __init__(self, name:str, desc:str) -> None:
        self.name = name
        self.desc = desc

    def set_value(self, value:float, unit:str):
        self.value = value
        self.unit = unit

    def set_comment(self, comment:str):
        self.comment = comment

    def use(self):
        print('{0} = {1} {2} ({3})'.format(self.name, self.value, self.unit, self.comment))


class ParamFormula:
    name:str
    desc:str
    value:float
    unit:str
    comment:str
    formula:str
    params:dict

    def __init__(self, name:str, desc:str) -> None:
        self.name = name
        self.desc = desc

    def set_unit(self, unit:str):
        self.unit = unit

    def set_comment(self, comment:str):
        self.comment = comment

    def set_formula(self, formula:str, params:dict):
        self.formula = formula
        self.params = params

    def use(self):
        self.value = ne.evaluate(self.formula, local_dict=self.params)
        self.comment = 'Вычисляем по формуле {0}'.format(self.formula)
        print('{0} = {1} {2} ({3})'.format(self.name, self.value, self.unit, self.comment))


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

        self.condenser = Condensator(self.compressor, self.agent, self.delta_t_v, self.t_in_v, self.t_out_v, self.t_k)

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

class Point:
    P:float = 0
    T:float = 0
    H:float = 0
    S:float = 0
    V:float = 0

    def set_p(self, p:float):
        self.P = p

    def set_t(self, t:float):
        self.T = t

    def set_h(self, h:float):
        self.H = h

    def set_s(self, s:float):
        self.S = s

    def set_v(self, v:float):
        self.V = v

    def get_params(self):
        return (self.P, self.T, self.H, self.S, self.V)
    
    def print_params(self):
        print(self.P, self.T, self.H, self.S, self.V)



class Graphic:

    solver:MainSolver
    dot_1_:Point
    dot_1:Point
    dot_2:Point
    dot_2_:Point
    dot_3_:Point
    dot_3:Point
    dot_4:Point

    fig = None

    def __init__(self, t_0:float, t_vs:float, t_k:float, t_p:float) -> None:

        self.t_0 = t_0
        self.t_vs = t_vs
        self.t_k = t_k
        self.t_p = t_p

        self.get_params_dot_1_()
        self.get_params_dot_1()
        self.get_params_dot_2()
        self.get_params_dot_2_()
        self.get_params_dot_3_()
        self.get_params_dot_3()
        self.get_params_dot_4()

    def make_graphic_template(self):
        p = []
        h = []

        conn = sqlite3.connect('main.db')
        c = conn.cursor()

        fig = go.Figure()

        # генерация линий x=0 и x=1
        p_l = []
        h_l = []
        hh_l = []

        c.execute('SELECT P,hd,hdd FROM amm_wet WHERE P<20 ORDER BY P ASC')

        res = c.fetchall()

        for r in res:
            p_l.append(r[0] / 10)
            h_l.append(r[1])
            hh_l.append(r[2])

        p_l.sort()
        h_l.sort()
        hh_l.sort()

        fig.add_trace(go.Scatter(x=h_l, y=p_l, name='x=0',
                            mode='lines', line=dict(color="black")))
        fig.add_trace(go.Scatter(x=hh_l, y=p_l, name='x=1',
                            mode='lines', line=dict(color="black")))

        # генерация изотерм на графике
        for x in range(203, 600, 5):
            p_l = []
            h_l = []

            c.execute('SELECT P,h FROM amm_water WHERE T={0} AND P<20 ORDER BY P DESC'.format(x))

            res = c.fetchall()

            for r in res:
                p_l.append(r[0] / 10)
                h_l.append(r[1])

            if x >= 200 and x <=400:
                c.execute('SELECT P,hd,hdd FROM amm_wet WHERE T={0} AND P<20 ORDER BY P DESC'.format(x))

                res = c.fetchall()

                for r in res:
                    p_l.append(r[0] / 10)
                    h_l.append(r[1])

                    p_l.append(r[0] / 10)
                    h_l.append(r[2])

            c.execute('SELECT P,h FROM amm_over WHERE T={0} AND P<20 ORDER BY P DESC'.format(x))

            res = c.fetchall()

            for r in res:
                p_l.append(r[0] / 10)
                h_l.append(r[1])
            
            fig.add_trace(go.Scatter(x=h_l, y=p_l, name='T={0}'.format(x-273), mode='lines+markers', line=dict(color="blue", dash="dash"), marker=dict(size=2)))

        # генерация изоэнтроп на графике
        for x in frange(8.8, 10.4, 0.1):
            mass = []
            data_type = [('P', float), ('h', float)]

            c.execute('SELECT P,h FROM amm_over WHERE s<' + str(x+0.009) + ' AND s>' + str(x-0.009) + ' AND P<20 ORDER BY P ASC')

            res = c.fetchall()

            for r in res:
                mass.append((float(r[0] / 10), float(r[1])))

            vdv = np.array(mass, dtype = data_type)

            vdv = np.sort(vdv, order = 'h')

            for n in range(0, len(vdv['h'])-1):
                if vdv['h'][n] >= vdv['h'][n+1]:
                    vdv['h'][n+1] = vdv['h'][n] + 0.0005

            z = np.polyfit(vdv['h'], vdv['P'], 5)
            f = np.poly1d(z)

            x_new = np.linspace(min(vdv['h']), max(vdv['h']), len(vdv['h']))
            y_new = f(x_new)

            fig.add_trace(go.Scatter(x=x_new, y=y_new, name='S={0}'.format(x),
                            mode='lines', line=dict(color="red")))

        # генерация изохор на графике                    
        for x in frange(0.1, 4, 0.1):
            mass = []
            data_type = [('P', float), ('h', float)]

            c.execute('SELECT P,h FROM amm_over WHERE v<' + str(x+0.007) + ' AND v>' + str(x-0.007) + ' AND P<20 ORDER BY P ASC')

            res = c.fetchall()

            for r in res:
                mass.append((float(r[0] / 10), float(r[1])))

            vdv = np.array(mass, dtype = data_type)

            vdv = np.sort(vdv, order = 'h')

            z = np.polyfit(vdv['h'], vdv['P'], 3)
            f = np.poly1d(z)

            x_new = np.linspace(min(vdv['h']), max(vdv['h']), len(vdv['h']))
            y_new = f(x_new)

            fig.add_trace(go.Scatter(x=x_new, y=y_new, name='v={0}'.format(x),
                                    mode='lines', line=dict(color="green")))
        
        fig.update_yaxes(type="log", ticklabelstep=2, showgrid=True, gridwidth=1, gridcolor='LightPink', zeroline=False, dtick = np.log10(1))
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightPink', zeroline=False, dtick = 50)

        self.fig = fig

    def make_graphic(self):
        x = [self.dot_1_.H, self.dot_1.H]
        y = [self.dot_1_.P, self.dot_1.P]
        self.fig.add_trace(go.Scatter(x=x, y=y, name='Перегрев пара перед компрессором', mode='lines+markers', line=dict(color='orange', width=5), marker=dict(color='orange', size=4)))

        x = [self.dot_1.H, self.dot_2.H]
        y = [self.dot_1.P, self.dot_2.P]
        self.fig.add_trace(go.Scatter(x=x, y=y, name='Сжатие в компрессоре', mode='lines+markers', line=dict(color='orange', width=5), marker=dict(color='orange', size=4)))

        x = [self.dot_2.H, self.dot_3.H]
        y = [self.dot_2.P, self.dot_3.P]
        self.fig.add_trace(go.Scatter(x=x, y=y, name='Охлаждение и конденсация в конденсаторе', mode='lines+markers', line=dict(color='orange', width=5), marker=dict(color='orange', size=4)))

        x = [self.dot_3.H, self.dot_4.H]
        y = [self.dot_3.P, self.dot_4.P]
        self.fig.add_trace(go.Scatter(x=x, y=y, name='Дросселирование', mode='lines+markers', line=dict(color='orange', width=5), marker=dict(color='orange', size=4)))

        x = [self.dot_4.H, self.dot_1_.H]
        y = [self.dot_4.P, self.dot_1_.P]
        self.fig.add_trace(go.Scatter(x=x, y=y, name='Кипение в испарителе', mode='lines+markers', line=dict(color='orange', width=5), marker=dict(color='orange', size=4)))

        self.fig.show()

    def get_params_dot_1_(self):
        temp_isp = self.t_0
        temp_isp_kelvin = temp_isp + 273

        c.execute('SELECT P,hdd,sdd FROM amm_wet WHERE T={0}'.format(round(temp_isp_kelvin)))

        res = c.fetchall()

        for r in res:
            self.dot_1_ = Point()
            self.dot_1_.set_p(r[0]/10)
            self.dot_1_.set_t(temp_isp)
            self.dot_1_.set_h(r[1])
            self.dot_1_.set_s(r[2])
            self.dot_1_.print_params()

    def get_params_dot_1(self):
        pressure_isp = self.dot_1_.P
        temp_insert = round(float(self.t_vs))
        temp_insert_kelvin = temp_insert + 273

        p = []
        c.execute('SELECT P FROM amm_over')
        res = c.fetchall()
        for r in res:
            p.append(float(r[0]))
        res = list(set(p))

        res.sort()

        pressure_val_min = min(res, key=lambda x:abs(x-pressure_isp*10))

        if pressure_val_min > pressure_isp:
            pressure_val_min = res[res.index(pressure_val_min) - 1]

        pressure_val_max = res[res.index(pressure_val_min) + 1]

        h_min = 0
        h_max = 0

        s_min = 0
        s_max = 0

        v_min = 0
        v_max = 0

        c.execute('SELECT h,s,v FROM amm_over WHERE P={0} AND T={1}'.format(pressure_val_min, temp_insert_kelvin))
        res = c.fetchall()
        for r in res:
            h_min = r[0]
            s_min = r[1]
            v_min = r[2]

        c.execute('SELECT h,s,v FROM amm_over WHERE P={0} AND T={1}'.format(pressure_val_max, temp_insert_kelvin))
        res = c.fetchall()
        for r in res:
            h_max = r[0]
            s_max = r[1]
            v_max = r[2]

        coeff = (pressure_isp*10-pressure_val_min)/(pressure_val_max-pressure_val_min)

        h_real = h_min + (h_max - h_min)*coeff
        s_real = s_min + (s_max - s_min)*coeff
        v_real = v_min + (v_max - v_min)*coeff

        self.dot_1 = Point()
        self.dot_1.set_p(pressure_isp)
        self.dot_1.set_t(float(self.t_vs))
        self.dot_1.set_h(h_real)
        self.dot_1.set_s(s_real)
        self.dot_1.set_v(v_real)

        self.dot_1.print_params()

    def get_params_dot_2(self):
        pressure_cond = 0
        pressure_cond_round = 0
        temp_cond = self.t_k
        temp_cond_kelvin = temp_cond + 273

        print(temp_cond_kelvin)

        c.execute('SELECT P FROM amm_wet WHERE T={0}'.format(round(temp_cond_kelvin)))

        res = c.fetchall()

        for r in res:
            pressure_cond = r[0]
            pressure_cond_round = round(r[0])

        c.execute('SELECT T,h,s FROM amm_over WHERE P={0} AND s>{1} AND s<{2}'.format(pressure_cond_round, self.dot_1.S - 0.005, self.dot_1.S + 0.005))

        res = c.fetchall()

        for r in res:
            self.dot_2 = Point()
            self.dot_2.set_p(pressure_cond / 10)
            self.dot_2.set_t(float(r[0]-273))
            self.dot_2.set_h(r[1])
            self.dot_2.set_s(r[2])
            self.dot_2.print_params()

    def get_params_dot_2_(self):
        temp_cond = self.t_k
        temp_cond_kelvin = temp_cond + 273

        c.execute('SELECT hdd,sdd FROM amm_wet WHERE T={0}'.format(round(temp_cond_kelvin)))

        res = c.fetchall()

        for r in res:
            self.dot_2_ = Point()
            self.dot_2_.set_p(self.dot_2.P)
            self.dot_2_.set_t(temp_cond)
            self.dot_2_.set_h(r[0])
            self.dot_2_.set_s(r[1])
            self.dot_2_.print_params()

    def get_params_dot_3_(self):
        temp_cond = self.t_k
        temp_cond_kelvin = temp_cond + 273

        c.execute('SELECT hd,sd FROM amm_wet WHERE T={0}'.format(round(temp_cond_kelvin)))

        res = c.fetchall()

        for r in res:
            self.dot_3_ = Point()
            self.dot_3_.set_p(self.dot_2.P)
            self.dot_3_.set_t(temp_cond)
            self.dot_3_.set_h(r[0])
            self.dot_3_.set_s(r[1])
            self.dot_3_.print_params()

    def get_params_dot_3(self):
        pressure_cond = self.dot_2.P
        temp_overcold = round(int(self.t_p))
        temp_overcold_kelvin = temp_overcold + 273

        p = []
        c.execute('SELECT P FROM amm_water')
        res = c.fetchall()
        for r in res:
            p.append(float(r[0]))
        res = list(set(p))

        res.sort()

        pressure_val_min = min(res, key=lambda x:abs(x-pressure_cond*10))

        if pressure_val_min > pressure_cond:
            pressure_val_min = res[res.index(pressure_val_min) - 1]

        pressure_val_max = res[res.index(pressure_val_min) + 1]

        if pressure_cond*10 == pressure_val_min:
            c.execute('SELECT h,s FROM amm_water WHERE P={0} AND T={1}'.format(pressure_val_min, temp_overcold_kelvin))
            res = c.fetchall()
            for r in res:
                h_min = r[0]
                s_min = r[1]

                self.dot_3 = Point()
                self.dot_3.set_p(pressure_cond)
                self.dot_3.set_t(self.t_p)
                self.dot_3.set_h(h_min)
                self.dot_3.set_s(s_min)

                self.dot_3.print_params()
        elif pressure_cond*10 == pressure_val_max:
            c.execute('SELECT h,s FROM amm_water WHERE P={0} AND T={1}'.format(pressure_val_max, temp_overcold_kelvin))
            res = c.fetchall()
            for r in res:
                h_max = r[0]
                s_max = r[1]

                self.dot_3 = Point()
                self.dot_3.set_p(pressure_cond)
                self.dot_3.set_t(self.t_p)
                self.dot_3.set_h(h_max)
                self.dot_3.set_s(s_max)

                self.dot_3.print_params()
        elif pressure_cond*10 < pressure_val_max and pressure_cond*10 > pressure_val_min:

            print(pressure_val_min, pressure_val_max)

            h_min = 0
            h_max = 0

            s_min = 0
            s_max = 0

            c.execute('SELECT h,s FROM amm_water WHERE P={0} AND T={1}'.format(pressure_val_min, temp_overcold_kelvin))
            res = c.fetchall()
            for r in res:
                h_min = r[0]
                s_min = r[1]

            c.execute('SELECT h,s FROM amm_water WHERE P={0} AND T={1}'.format(pressure_val_max, temp_overcold_kelvin))
            res = c.fetchall()
            for r in res:
                h_max = r[0]
                s_max = r[1]

            coeff = (pressure_cond*10-pressure_val_min)/(pressure_val_max-pressure_val_min)

            h_real = h_min + (h_max - h_min)*coeff
            s_real = s_min + (s_max - s_min)*coeff

            self.dot_3 = Point()
            self.dot_3.set_p(pressure_cond)
            self.dot_3.set_t(self.t_p)
            self.dot_3.set_h(h_real)
            self.dot_3.set_s(s_real)

            self.dot_3.print_params()

    def get_params_dot_4(self):
        self.dot_4 = Point()
        self.dot_4.set_p(self.dot_1_.P)
        self.dot_4.set_t(self.dot_1_.T)
        self.dot_4.set_h(self.dot_3.H)
        self.dot_4.set_s(0)
        self.dot_4.print_params()

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

class Compressor:

    mark:str
    diameter:float
    piston_stroke:float
    theor_v: float
    frequency:int
    Q_R717:float
    N_R717:float

    def __init__(self, params:dict) -> None:
        self.mark = params['mark']
        self.diameter = params['diameter']
        self.piston_stroke = params['piston_stroke']
        self.theor_v = params['theor_v']
        self.frequency = params['frequency']
        self.Q_R717 = params['Q_R717']
        self.N_R717 = params['N_R717']

class Condensator:
    
    compressor:Compressor
    coolant:str
    K:Param
    Qk:ParamFormula
    Cp:Param
    tv_in:Param
    tv_out:ParamFormula
    tk:ParamFormula
    tv_delta:Param

    def __init__(self, compressor:Compressor, coolant:str, delta_t_v:Param, t_in_v:Param, t_out_v:ParamFormula, t_k:ParamFormula) -> None:
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
        self.Gv.set_formula('Qk/(Cp*delta_tv)', {'Qk': self.Qk.value*1000, 'Cp': self.Cp.value, 'delta_tv': self.delta_t_v.value})
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
        self.delta_t_log.set_formula('tmax - tmin/()', {})
        self.delta_t_log.set_unit('кг/с')
        self.delta_t_log.use()

        


        

solver = MainSolver(180.0, 'Владимир', 'R717', -2.0, 'КТ', 'оборотное')

    
