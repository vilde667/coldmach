import plotly.graph_objs as go
import numpy as np
from point import Point
from funcs import frange

from tables import R717Liquid, R717Wet, R717Over


class Graphic:
    dot_1_: Point
    dot_1: Point
    dot_2: Point
    dot_2_: Point
    dot_3_: Point
    dot_3: Point
    dot_4: Point
    show_graphic: bool

    fig = None

    def __init__(self, t_0: float, t_vs: float, t_k: float, t_p: float, show_graphic: bool) -> None:

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

        self.show_graphic = show_graphic

    def make_graphic_template(self):
        p = []
        h = []

        fig = go.Figure()

        # генерация линий x=0 и x=1
        p_l = []
        h_l = []
        hh_l = []

        results = R717Wet.select().where(R717Wet.P < 20).order_by(R717Wet.P.asc())
        for r in results:
            p_l.append(r.P / 10)
            h_l.append(r.hd)
            hh_l.append(r.hdd)

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

            results = R717Liquid.select().where(R717Liquid.T == x).where(R717Liquid.P < 20).order_by(R717Liquid.P.desc())
            for r in results:
                p_l.append(r.P / 10)
                h_l.append(r.h)

            if 200 <= x <= 400:
                results = R717Wet.select().where(R717Wet.T == x).where(R717Wet.P < 20).order_by(R717Wet.P.desc())
                for r in results:
                    p_l.append(r.P / 10)
                    h_l.append(r.hd)

                    p_l.append(r.P / 10)
                    h_l.append(r.hdd)

            results = R717Over.select().where(R717Over.T == x).where(R717Over.P < 20).order_by(R717Over.P.desc())
            for r in results:
                p_l.append(r.P / 10)
                h_l.append(r.h)

            fig.add_trace(go.Scatter(x=h_l, y=p_l, name='T={0}'.format(x - 273), mode='lines+markers',
                                     line=dict(color="blue", dash="dash"), marker=dict(size=2)))

        # генерация изоэнтроп на графике
        for x in frange(8.8, 10.4, 0.1):
            mass = []
            data_type = [('P', float), ('h', float)]

            results = R717Over.select().where(R717Over.s < x + 0.009).where(R717Over.s > x - 0.009)\
                .where(R717Over.P < 20).order_by(R717Over.P.asc())
            for r in results:
                mass.append((float(r.P / 10), float(r.h)))

            vdv = np.array(mass, dtype=data_type)

            vdv = np.sort(vdv, order='h')

            for n in range(0, len(vdv['h']) - 1):
                if vdv['h'][n] >= vdv['h'][n + 1]:
                    vdv['h'][n + 1] = vdv['h'][n] + 0.0005

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

            results = R717Over.select().where(R717Over.v < x + 0.007).where(R717Over.v > x - 0.007) \
                .where(R717Over.P < 20).order_by(R717Over.P.asc())
            for r in results:
                mass.append((float(r.P / 10), float(r.h)))

            vdv = np.array(mass, dtype=data_type)

            vdv = np.sort(vdv, order='h')

            z = np.polyfit(vdv['h'], vdv['P'], 3)
            f = np.poly1d(z)

            x_new = np.linspace(min(vdv['h']), max(vdv['h']), len(vdv['h']))
            y_new = f(x_new)

            fig.add_trace(go.Scatter(x=x_new, y=y_new, name='v={0}'.format(x),
                                     mode='lines', line=dict(color="green")))

        fig.update_yaxes(type="log", ticklabelstep=2, showgrid=True, gridwidth=1, gridcolor='LightPink', zeroline=False,
                         dtick=np.log10(1))
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightPink', zeroline=False, dtick=50)

        self.fig = fig

    def make_graphic(self):
        x = [self.dot_1_.H, self.dot_1.H]
        y = [self.dot_1_.P, self.dot_1.P]
        self.fig.add_trace(go.Scatter(x=x, y=y, name='Перегрев пара перед компрессором', mode='lines+markers',
                                      line=dict(color='orange', width=5), marker=dict(color='orange', size=4)))

        x = [self.dot_1.H, self.dot_2.H]
        y = [self.dot_1.P, self.dot_2.P]
        self.fig.add_trace(
            go.Scatter(x=x, y=y, name='Сжатие в компрессоре', mode='lines+markers', line=dict(color='orange', width=5),
                       marker=dict(color='orange', size=4)))

        x = [self.dot_2.H, self.dot_3.H]
        y = [self.dot_2.P, self.dot_3.P]
        self.fig.add_trace(go.Scatter(x=x, y=y, name='Охлаждение и конденсация в конденсаторе', mode='lines+markers',
                                      line=dict(color='orange', width=5), marker=dict(color='orange', size=4)))

        x = [self.dot_3.H, self.dot_4.H]
        y = [self.dot_3.P, self.dot_4.P]
        self.fig.add_trace(
            go.Scatter(x=x, y=y, name='Дросселирование', mode='lines+markers', line=dict(color='orange', width=5),
                       marker=dict(color='orange', size=4)))

        x = [self.dot_4.H, self.dot_1_.H]
        y = [self.dot_4.P, self.dot_1_.P]
        self.fig.add_trace(
            go.Scatter(x=x, y=y, name='Кипение в испарителе', mode='lines+markers', line=dict(color='orange', width=5),
                       marker=dict(color='orange', size=4)))

        if self.show_graphic:
            self.fig.show()

    def get_params_dot_1_(self):
        temp_isp = self.t_0
        temp_isp_kelvin = temp_isp + 273

        result = R717Wet.select().where(R717Wet.T == round(temp_isp_kelvin)).get()
        self.dot_1_ = Point()
        self.dot_1_.set_p(result.P / 10)
        self.dot_1_.set_t(temp_isp)
        self.dot_1_.set_h(result.hdd)
        self.dot_1_.set_s(result.sdd)
        self.dot_1_.print_params()

    def get_params_dot_1(self):
        pressure_isp = self.dot_1_.P
        temp_insert = round(float(self.t_vs))
        temp_insert_kelvin = temp_insert + 273

        p = []
        results = R717Over.select()
        for r in results:
            p.append(float(r.P))
        res = list(set(p))

        res.sort()

        pressure_val_min = min(res, key=lambda x: abs(x - pressure_isp * 10))

        if pressure_val_min > pressure_isp:
            pressure_val_min = res[res.index(pressure_val_min) - 1]

        pressure_val_max = res[res.index(pressure_val_min) + 1]

        h_min = 0
        h_max = 0

        s_min = 0
        s_max = 0

        v_min = 0
        v_max = 0

        results = R717Over.select().where(R717Over.T == temp_insert_kelvin).where(R717Over.P == pressure_val_min)
        for r in results:
            h_min = r.h
            s_min = r.s
            v_min = r.v

        results = R717Over.select().where(R717Over.T == temp_insert_kelvin).where(R717Over.P == pressure_val_max)
        for r in results:
            h_max = r.h
            s_max = r.s
            v_max = r.v

        coeff = (pressure_isp * 10 - pressure_val_min) / (pressure_val_max - pressure_val_min)

        h_real = h_min + (h_max - h_min) * coeff
        s_real = s_min + (s_max - s_min) * coeff
        v_real = v_min + (v_max - v_min) * coeff

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

        results = R717Wet.select().where(R717Wet.T == round(temp_cond_kelvin))

        for r in results:
            pressure_cond = r.P
            pressure_cond_round = round(r.P)

        results = R717Over.select().where(R717Over.P == pressure_cond_round)\
            .where(R717Over.s > self.dot_1.S - 0.005).where(R717Over.s < self.dot_1.S + 0.005)

        for r in results:
            self.dot_2 = Point()
            self.dot_2.set_p(pressure_cond / 10)
            self.dot_2.set_t(float(r.T - 273))
            self.dot_2.set_h(r.h)
            self.dot_2.set_s(r.s)
            self.dot_2.print_params()

    def get_params_dot_2_(self):
        temp_cond = self.t_k
        temp_cond_kelvin = temp_cond + 273

        results = R717Wet.select().where(R717Wet.T == round(temp_cond_kelvin))

        for r in results:
            self.dot_2_ = Point()
            self.dot_2_.set_p(self.dot_2.P)
            self.dot_2_.set_t(temp_cond)
            self.dot_2_.set_h(r.hdd)
            self.dot_2_.set_s(r.sdd)
            self.dot_2_.print_params()

    def get_params_dot_3_(self):
        temp_cond = self.t_k
        temp_cond_kelvin = temp_cond + 273

        results = R717Wet.select().where(R717Wet.T == round(temp_cond_kelvin))

        for r in results:
            self.dot_3_ = Point()
            self.dot_3_.set_p(self.dot_2.P)
            self.dot_3_.set_t(temp_cond)
            self.dot_3_.set_h(r.hd)
            self.dot_3_.set_s(r.sd)
            self.dot_3_.print_params()

    def get_params_dot_3(self):
        pressure_cond = self.dot_2.P
        temp_overcold = round(int(self.t_p))
        temp_overcold_kelvin = temp_overcold + 273

        p = []
        results = R717Liquid.select()
        for r in results:
            p.append(float(r.P))
        res = list(set(p))

        res.sort()

        pressure_val_min = min(res, key=lambda x: abs(x - pressure_cond * 10))

        if pressure_val_min > pressure_cond:
            pressure_val_min = res[res.index(pressure_val_min) - 1]

        pressure_val_max = res[res.index(pressure_val_min) + 1]

        if pressure_cond * 10 == pressure_val_min:
            results = R717Liquid.select().where(R717Liquid.P == pressure_val_min)\
                .where(R717Liquid.T == temp_overcold_kelvin)
            for r in results:
                h_min = r.h
                s_min = r.s

                self.dot_3 = Point()
                self.dot_3.set_p(pressure_cond)
                self.dot_3.set_t(self.t_p)
                self.dot_3.set_h(h_min)
                self.dot_3.set_s(s_min)

                self.dot_3.print_params()
        elif pressure_cond * 10 == pressure_val_max:
            results = R717Liquid.select().where(R717Liquid.P == pressure_val_max) \
                .where(R717Liquid.T == temp_overcold_kelvin)
            for r in results:
                h_max = r.h
                s_max = r.s

                self.dot_3 = Point()
                self.dot_3.set_p(pressure_cond)
                self.dot_3.set_t(self.t_p)
                self.dot_3.set_h(h_max)
                self.dot_3.set_s(s_max)

                self.dot_3.print_params()
        elif pressure_val_max > pressure_cond * 10 > pressure_val_min:

            print(pressure_val_min, pressure_val_max)

            h_min = 0
            h_max = 0

            s_min = 0
            s_max = 0

            results = R717Liquid.select().where(R717Liquid.P == pressure_val_min) \
                .where(R717Liquid.T == temp_overcold_kelvin)
            for r in results:
                h_min = r.h
                s_min = r.s

            results = R717Liquid.select().where(R717Liquid.P == pressure_val_max) \
                .where(R717Liquid.T == temp_overcold_kelvin)
            for r in results:
                h_max = r.h
                s_max = r.s

            coeff = (pressure_cond * 10 - pressure_val_min) / (pressure_val_max - pressure_val_min)

            h_real = h_min + (h_max - h_min) * coeff
            s_real = s_min + (s_max - s_min) * coeff

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
