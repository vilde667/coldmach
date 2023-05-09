import numpy as np
import math
import sqlite3

from params import Param
from params import ParamFormula

from compressor.compressor import Compressor
from condenser.condenser import Condenser

from graphic import Graphic

from tables import Cities, Compressors, FeedRatio, get_approx_value


class MainSolver:
    compressor: Compressor
    delta_t: Param
    delta_t_p: Param
    t_brine_in: ParamFormula
    delta_t_v: Param
    t_in_v: Param
    t_out_v: ParamFormula
    t_0: ParamFormula
    delta_t_k: Param
    t_k: ParamFormula
    show_graphics: bool

    def __init__(self, Q0, city, agent, temp_brine_out, condenser_type, water_supply_type, show_graphics: bool) -> None:
        self.Q0 = Param('Q0', 'Расчетная холодопроизводительность')
        self.Q0.set_value(Q0, 'кВт')

        self.city = city

        self.agent = agent

        self.t_brine_out = Param('t\'\'р', 'Температура рассола на выходе')
        self.t_brine_out.set_value(temp_brine_out, '°С')

        self.condenser_type = condenser_type

        self.water_supply_type = water_supply_type

        self.show_graphics = show_graphics

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

        self.graphic = Graphic(self.t_0.value, self.t_vs.value, self.t_k.value, self.t_p.value, self.show_graphics)

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

        self.standart_graphic = Graphic(-15, -10, 30, 25, self.show_graphics)

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

        res = Cities.select().where(Cities.city == city).get()
        result = res.t_air_out

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
        value = get_approx_value(FeedRatio, FeedRatio.R717, 'R717', FeedRatio.compression_ratio, 'compression_ratio',
                                 self.compression_ratio.value)
        self.lambda_real = Param('λ', 'коэффициент подачи реального компрессора')
        self.lambda_real.set_value(value, '')
        self.lambda_real.set_comment('Находим по графику 4.11 для R717 по степени сжатия={0}'.format(self.compression_ratio.value))
        self.lambda_real.use()

    def set_lambda_standart(self):
        value = get_approx_value(FeedRatio, FeedRatio.R717, 'R717', FeedRatio.compression_ratio, 'compression_ratio',
                                 self.compression_ratio_standart.value)
        self.lambda_standart = Param('λ', 'коэффициент подачи стандартного компрессора')
        self.lambda_standart.set_value(value, '')
        self.lambda_standart.set_comment('Находим по графику 4.11 для R717 по степени сжатия={0}'.format(self.compression_ratio_standart.value))
        self.lambda_standart.use()

    def set_Q0_st(self):
        self.Q0_st = ParamFormula('Q0 ст', 'стандартная холодопроизводительность компрессора')
        self.Q0_st.set_formula('Q0 * (q_v_st * lambda_standart)/(q_v * lambda_real)', {'Q0': self.Q0.value, 'q_v_st': self.q_v_st.value, 'q_v': self.q_v.value, 'lambda_standart': self.lambda_standart.value, 'lambda_real': self.lambda_real.value})
        self.Q0_st.set_unit('кВт')
        self.Q0_st.use()

    def set_compressor(self):

        compressor = Compressors.select().where(Compressors.Q_R717 > self.Q0_st.value).order_by(Compressors.Q_R717)\
                     .get()

        params = {'mark': compressor.mark, 'diameter': compressor.diameter, 'piston_stroke': compressor.piston_stroke,
                  'theor_v': compressor.theor_v, 'frequency': compressor.frequency, 'Q_R717': compressor.Q_R717,
                  'N_R717': compressor.N_R717}

        self.compressor = Compressor(params)

        print('По полученной холодопроизводительности выбираем компрессор марки ' + compressor.mark)
    
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


solver = MainSolver(180.0, 'Владимир', 'R717', -2.0, 'КТ', 'оборотное', True)
