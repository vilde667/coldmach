import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np

from point import Point

from funcs import frange

class Graphic:

    solver: MainSolver
    dot_1_: Point
    dot_1: Point
    dot_2: Point
    dot_2_: Point
    dot_3_: Point
    dot_3: Point
    dot_4: Point

    fig = None

    def __init__(self, t_0: float, t_vs: float, t_k: float, t_p: float) -> None:

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