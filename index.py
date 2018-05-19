import sqlite3
import matplotlib.pyplot as plt
import numpy as np


class Value:
    def __init__(self, unit, value, amount):
        self.amount = amount
        self.value = float(value)
        self.unit = unit


class Point:
    def __init__(self):
        self.temp = Value('градус', 0, 'температура')
        self.press = Value('МПа', 0, 'давление')
        self.ent = Value('кДж/кг', 0, 'энтальпия')
        self.cap = Value('м3/кг', 0, 'удельный объем')
        self.state = ''
        self.human_name = ''

    def print_all(self):
        print('Для данной точки', self.human_name)
        for y in [self.temp, self.press, self.ent, self.cap]:
            self.get_one(y)

    @staticmethod
    def get_one(param):
        print('{0} = {1} {2}'.format(param.amount, param.value, param.unit))

    @staticmethod
    def get_state(obj):
        r = input('''введите состояние хлодоагента в данной точке\n
                 wet - влажный насыщенный пар\n
                 over - перегретый пар\n
                 liquid - жидкость''')
        if r == 'liquid':
            obj.cap.value = 0
        return r


class StandartParams:
    def __init__(self):
        self.t_0 = -15
        self.t_vs = -10
        self.t_k = 30
        self.t_p = 25
        self.point1 = Point()
        self.point1d = Point()
        self.point2 = Point()
        self.point2d = Point()
        self.point3 = Point()
        self.point3d = Point()
        self.point4 = Point()
        self.point4d = Point()
        self.points = []


class ColdSystem:
    def __init__(self, q0, ca, city, t_out_r):
        """q0 - холодопроизводительность(int), Ca - хладоагент(str),
           city - город установки(str), t_out_r - температура рассола на выходе из испарителя(float)"""
        self.Q0 = q0
        self.Ca = ca
        self.city = city
        self.t = {'out_r': t_out_r}
        self.point1 = Point()
        self.point1d = Point()
        self.point2 = Point()
        self.point2d = Point()
        self.point3 = Point()
        self.point3d = Point()
        self.point4 = Point()
        self.point4d = Point()
        self.std = StandartParams()

    @staticmethod
    def __get_t_in_k(city):
        conn = sqlite3.connect('main.db')
        c = conn.cursor()
        t = (city,)
        c.execute('SELECT * FROM citys WHERE city=?', t)
        response = c.fetchone()
        return float(response[1].replace(',', '.'))

    def __get_points_4_4d(self, t, h_4, h_4d):
        hd, hdd, vd, vdd = self.__query('amm_wet', 'hd,hdd,vd,vdd', t=t)
        x_4 = (h_4 - hd) / (hdd - hd)
        v_4 = vd + x_4 * (vdd - vd)
        x_4d = (h_4d - hd) / (hdd - hd)
        v_4d = vd + x_4d * (vdd - vd)
        return v_4, v_4d

    def __get_point3(self, p=None, t=None):
        t_tab_min, t_tab_max = self.__get_near_t(t)
        p_tab_min, p_tab_max = self.__get_near_p(p)
        res_all = [self.__query('amm_over', 'v,h', p=p_tab_min, t=int(t_tab_min)),
                   self.__query('amm_over', 'v,h', p=p_tab_max, t=int(t_tab_min))]
        zipped = tuple(zip(res_all[0], res_all[1]))
        v = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[0][1] - zipped[0][0]) + zipped[0][0]
        h = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[1][1] - zipped[1][0]) + zipped[1][0]
        return h, v

    def __get_wet_p_max(self, p=None):
        v, h = self.__query('amm_wet', 'vd,hd', p=p)
        return v, h

    def __get_wet_p_min(self, p=None):
        v, h = self.__query('amm_wet', 'vd,hd', p=p)
        return v, h

    def __get_near_p(self, p=None):
        delta_p = 10
        response = self.__query('amm_over', 'P', p1=p - delta_p, p2=p + delta_p)
        near_minus = []
        near_plus = []
        for x in response:
            if p - x[0] >= 0:
                near_plus.append(p - x[0])
            else:
                near_minus.append(p - x[0])
        p_tab_min_ = max(near_minus)
        p_tab_max_ = min(near_plus)
        p_tab_min = p - p_tab_max_
        p_tab_max = p - p_tab_min_
        return p_tab_min, p_tab_max

    def __get_near_t(self, t=None):
        delta_t = 10
        response = self.__query('amm_over', 'T', t1=t - delta_t, t2=t + delta_t)
        near_minus = []
        near_plus = []
        for x in response:
            if t - x[0] >= 0:
                near_plus.append(t - x[0])
            else:
                near_minus.append(t - x[0])
        t_tab_min_ = max(near_minus)
        t_tab_max_ = min(near_plus)
        t_tab_min = t - t_tab_min_
        t_tab_max = t - t_tab_max_
        return t_tab_max, t_tab_min

    def __get_over_p_t(self, p=None, t=None, p_tab_min=None, p_tab_max=None):
        print(t, p)
        t_tab_min, t_tab_max = self.__get_near_t(t)

        res_all = [self.__query('amm_over', 'v,h,s', p=p_tab_min, t=t_tab_min),
                   self.__query('amm_over', 'v,h,s', p=p_tab_max, t=t_tab_min),
                   self.__query('amm_over', 'v,h,s', p=p_tab_min, t=t_tab_max),
                   self.__query('amm_over', 'v,h,s', p=p_tab_max, t=t_tab_max)]
        zipped = tuple(zip(res_all[0], res_all[1], res_all[2], res_all[3]))
        v1 = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[0][2] - zipped[0][0]) + zipped[0][0]
        v2 = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[0][3] - zipped[0][1]) + zipped[0][1]
        v = (p - p_tab_min) / (p_tab_max - p_tab_min) * (v2 - v1) + v1
        h1 = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[1][2] - zipped[1][0]) + zipped[1][0]
        h2 = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[1][3] - zipped[1][1]) + zipped[1][1]
        h = (p - p_tab_min) / (p_tab_max - p_tab_min) * (h2 - h1) + h1
        s1 = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[2][2] - zipped[2][0]) + zipped[2][0]
        s2 = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[2][3] - zipped[2][1]) + zipped[2][1]
        s = (p - p_tab_min) / (p_tab_max - p_tab_min) * (s2 - s1) + s1
        print(zipped)
        print(v, h, s)
        return v, h, s, p_tab_min, p_tab_max

    def __get_over_p_s(self, p=None, s=None, p_tab_min=None, p_tab_max=None):
        print(s, p)
        delta_s = 0.2
        response = self.__query('amm_over', 't,s', p=p_tab_min, s1=s - delta_s, s2=s + delta_s)
        near_plus = []
        near_min = []
        for x in response:
            if s - x[1] >= 0:
                near_plus.append(s - x[1])
            else:
                near_min.append(s - x[1])
        s_tab_min = max(near_min)
        s_tab_max = min(near_plus)
        s_tab_min1 = s - s_tab_min
        s_tab_max1 = s - s_tab_max
        response = self.__query('amm_over', 't,s', p=p_tab_max, s1=s - delta_s, s2=s + delta_s)
        near_plus = []
        near_min = []
        for x in response:
            if s - x[1] >= 0:
                near_plus.append(s - x[1])
            else:
                near_min.append(s - x[1])
        s_tab_min = max(near_min)
        s_tab_max = min(near_plus)
        s_tab_min2 = s - s_tab_min
        s_tab_max2 = s - s_tab_max
        print(s_tab_min1, s_tab_min2, s_tab_max1, s_tab_max2)
        res_all = [self.__query('amm_over', 't,h,v', p=p_tab_min, s=s_tab_min1),
                   self.__query('amm_over', 't,h,v', p=p_tab_min, s=s_tab_max1),
                   self.__query('amm_over', 't,h,v', p=p_tab_max, s=s_tab_min2),
                   self.__query('amm_over', 't,h,v', p=p_tab_max, s=s_tab_max2)]
        zipped = tuple(zip(res_all[0], res_all[1], res_all[2], res_all[3]))
        t1 = (s - s_tab_min1) / (s_tab_max1 - s_tab_min1) * (zipped[0][1] - zipped[0][0]) + zipped[0][0]
        t2 = (s - s_tab_min2) / (s_tab_max2 - s_tab_min2) * (zipped[0][3] - zipped[0][2]) + zipped[0][2]
        t = (p - p_tab_min) / (p_tab_max - p_tab_min) * (t2 - t1) + t1
        h1 = (s - s_tab_min1) / (s_tab_max1 - s_tab_min1) * (zipped[1][1] - zipped[1][0]) + zipped[1][0]
        h2 = (s - s_tab_min2) / (s_tab_max2 - s_tab_min2) * (zipped[1][3] - zipped[1][2]) + zipped[1][2]
        h = (p - p_tab_min) / (p_tab_max - p_tab_min) * (h2 - h1) + h1
        v1 = (s - s_tab_min1) / (s_tab_max1 - s_tab_min1) * (zipped[2][1] - zipped[2][0]) + zipped[2][0]
        v2 = (s - s_tab_min2) / (s_tab_max2 - s_tab_min2) * (zipped[2][3] - zipped[2][2]) + zipped[2][2]
        v = (p - p_tab_min) / (p_tab_max - p_tab_min) * (v2 - v1) + v1
        return t, h, v

    @staticmethod
    def __get_lambda(pk, p0):
        delta = pk / p0
        func = [(4.3, 0.8), (5, 0.79), (6, 0.77), (7, 0.755), (8, 0.745), (9, 0.73)]
        for x in range(5):
            if func[x][0] < delta < func[x + 1][0]:
                lambda_ = (delta - func[x][0]) / (func[x + 1][0] - func[x][0]) * (func[x + 1][1] - func[x][1]) + \
                          func[x][1]
                return lambda_

    @staticmethod
    def __query(table, sel, **kwargs):
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        try:
            if 'p' in kwargs.keys() and 's' in kwargs.keys():
                cur.execute('SELECT {0} FROM {1} WHERE P={2} AND S={3};'.format(sel, table, kwargs['p'], kwargs['s']))
                r = cur.fetchone()
                return r
            elif 'p' in kwargs.keys() and 's1' in kwargs.keys() and 's2' in kwargs.keys():
                cur.execute('SELECT {0} FROM {1} WHERE P={2} AND S>{3} AND S<{4};'.format(
                    sel, table, kwargs['p'], kwargs['s1'], kwargs['s2']))
                r = cur.fetchall()
                return r
            elif 'p' in kwargs.keys() and 't' in kwargs.keys():
                cur.execute('SELECT {0} FROM {1} WHERE P={2} AND T={3};'.format(
                    sel, table, kwargs['p'], kwargs['t']))
                r = cur.fetchone()
                return r
            elif 't1' in kwargs.keys() and 't2' in kwargs.keys():
                cur.execute('SELECT {0} FROM {1} WHERE P=3.0 AND T>{2} AND T<{3};'.format(
                    sel, table, kwargs['t1'], kwargs['t2']))
                r = cur.fetchall()
                return r
            elif 'p1' in kwargs.keys() and 'p2' in kwargs.keys():
                cur.execute('SELECT {0} FROM {1} WHERE T=300.0 AND P>{2} AND P<{3};'.format(
                    sel, table, kwargs['p1'], kwargs['p2']))
                r = cur.fetchall()
                return r
            elif 'p' in kwargs.keys() and 's' in kwargs.keys():
                cur.execute('SELECT {0} FROM {1} WHERE P={2} AND S={3};'.format(sel, table, kwargs['p'], kwargs['s']))
                r = cur.fetchone()
                return r
            elif 'p' in kwargs.keys():
                cur.execute('SELECT {0} FROM {1} WHERE P={2};'.format(sel, table, kwargs['p']))
                r = cur.fetchone()
                return r
            elif 't' in kwargs.keys():
                cur.execute('SELECT {0} FROM {1} WHERE T={2};'.format(sel, table, kwargs['t']))
                r = cur.fetchone()
                return r
            else:
                return 0
        finally:
            conn.close()

    def get_std_params(self):
        std = StandartParams()  # type: StandartParams
        # точка 1
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT * FROM amm_wet WHERE T={0};'.format(int(std.t_0 + 273.15)))
        req = cur.fetchone()
        p_1 = req[2]
        p_1_tab_min, p_1_tab_max = self.__get_near_p(p_1)

        v_1, h_1, s_1, p_tab_min, p_tab_max = self.__get_over_p_t(p_1, std.t_vs + 273.15, p_1_tab_min, p_1_tab_max)
        std.point1.temp.value = std.t_vs
        std.point1.press.value = p_1 / 10
        std.point1.ent.value = h_1
        std.point1.cap.value = v_1
        std.point1.state = 'over'
        std.point1.human_name = '1'
        std.point1.print_all()

        # точка 1'
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT vdd, hdd FROM amm_wet WHERE T={0};'.format(round(std.t_0 + 273.15)))
        res_1d = cur.fetchone()
        conn.close()
        v_1d, h_1d = res_1d
        std.point1d.temp.value = std.t_0
        std.point1d.press.value = p_1 / 10
        std.point1d.ent.value = h_1d
        std.point1d.cap.value = v_1d
        std.point1d.state = 'wet'
        std.point1d.human_name = '1\''
        std.point1d.print_all()

        # точка 2
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT * FROM amm_wet WHERE T={0};'.format(round(std.t_k + 273.15)))
        req = cur.fetchone()
        conn.close()
        p_2 = req[2]
        p_2_tab_min, p_2_tab_max = self.__get_near_p(p_2)
        t_2, h_2, v_2 = self.__get_over_p_s(p_2, s_1, p_2_tab_min, p_2_tab_max)

        std.point2.temp.value = t_2 - 273.15
        std.point2.press.value = p_2 / 10
        std.point2.ent.value = h_2
        std.point2.cap.value = v_2
        std.point2.state = 'over'
        std.point2.human_name = '2'
        std.point2.print_all()

        # точка 2'
        v_2d, h_2d = self.__get_wet_p_max(p_2)
        std.point2d.temp.value = std.t_k
        std.point2d.press.value = std.point2.press.value
        std.point2d.ent.value = h_2d
        std.point2d.cap.value = v_2d
        std.point2d.state = 'wet'
        std.point2d.human_name = '2\''
        std.point2d.print_all()

        # точка 3'
        v_3d, h_3d = self.__get_wet_p_min(p_2)
        std.point3d.temp.value = std.t_k
        std.point3d.press.value = std.point2.press.value
        std.point3d.ent.value = h_3d
        std.point3d.cap.value = v_3d
        std.point3d.state = 'wet'
        std.point3d.human_name = '3\''
        std.point3d.print_all()

        # точка 3
        h_3, v_3 = self.__get_point3(std.point2d.press.value * 10, std.t_p + 273.15)
        std.point3.temp.value = std.t_p
        std.point3.press.value = std.point2.press.value
        std.point3.ent.value = h_3
        std.point3.cap.value = v_3
        std.point3.state = 'liquid'
        std.point3.human_name = '3'
        std.point3.print_all()

        # точка 4
        v_4, v_4d = self.__get_points_4_4d(round(std.t_0 + 273.15), h_3, h_3d)

        std.point4.temp.value = std.t_0
        std.point4.press.value = std.point1.press.value
        std.point4.ent.value = std.point3.ent.value
        std.point4.cap.value = v_4
        std.point4.state = 'wet'
        std.point4.human_name = '4'
        std.point4.print_all()

        # точка 4'
        std.point4d.temp.value = std.t_0
        std.point4d.press.value = std.point1.press.value
        std.point4d.ent.value = std.point3d.ent.value
        std.point4d.cap.value = v_4d
        std.point4d.state = 'wet'
        std.point4d.human_name = '4\''
        std.point4d.print_all()

        x = np.arange(std.point1.ent.value, std.point2.ent.value)
        p = std.point2.press.value
        r = std.point2.ent.value
        g = (p - std.point1.press.value) / ((std.point1.ent.value - r) ** 2)
        z = -g * (x - r) ** 2 + p
        std.points = [std.point1, std.point1d, std.point2, std.point2d, std.point3, std.point3d, std.point4,
                      std.point4d]
        plt.plot([std.point2.ent.value, std.point3.ent.value], [std.point2.press.value, std.point3.press.value])
        plt.plot([std.point3.ent.value, std.point4.ent.value], [std.point3.press.value, std.point4.press.value])
        plt.plot([std.point4.ent.value, std.point1.ent.value], [std.point4.press.value, std.point1.press.value])
        for point in std.points:
            plt.scatter(point.ent.value, point.press.value)
            if '\'' in point.human_name:
                plt.annotate(s=point.human_name, xy=(point.ent.value, point.press.value,),
                             xytext=(point.ent.value, point.press.value + 0.03,))
            else:
                plt.annotate(s=point.human_name, xy=(point.ent.value, point.press.value,),
                             xytext=(point.ent.value + 10, point.press.value - 0.06,))
        plt.plot(x, z)
        plt.xlim(200, 2100)
        plt.ylim(0.1, 1.5)
        plt.show()

        self.std = std

    def plot_graph(self):
        x = np.arange(self.point1.ent.value, self.point2.ent.value)
        p = self.point2.press.value
        r = self.point2.ent.value
        g = (p - self.point1.press.value) / ((self.point1.ent.value - r) ** 2)
        z = -g * (x - r) ** 2 + p
        plt.plot([self.point2.ent.value, self.point3.ent.value], [self.point2.press.value, self.point3.press.value])
        plt.plot([self.point3.ent.value, self.point4.ent.value], [self.point3.press.value, self.point4.press.value])
        plt.plot([self.point4.ent.value, self.point1.ent.value], [self.point4.press.value, self.point1.press.value])
        for point in (self.point1, self.point1d, self.point2, self.point2d, self.point3, self.point3d, self.point4,
                      self.point4d):
            plt.scatter(point.ent.value, point.press.value)
            if '\'' in point.human_name:
                plt.annotate(s=point.human_name, xy=(point.ent.value, point.press.value,),
                             xytext=(point.ent.value, point.press.value + 0.03,))
            else:
                plt.annotate(s=point.human_name, xy=(point.ent.value, point.press.value,),
                             xytext=(point.ent.value + 10, point.press.value - 0.06,))
        plt.plot(x, z)
        plt.xlim(200, 2100)
        plt.ylim(0.1, 1.5)
        plt.show()

    def get_charact_params(self):
        pass
        # q0 = self.point1d.ent.value - self.point4.ent.value
        # l = self.point2.ent.value - self.point1.ent.value
        # qk = self.point2.ent.value - self.point3d.ent.value
        # e = q0 / l
        # g0 = self.Q0 / q0
        # v0 = g0 * self.point1.cap.value
        # qv = q0 / self.point1.cap.value
        # qv_st = (self.std.point1d.ent.value - self.std.point4.ent.value) / self.std.point1.cap.value
        # lambda_ = self.__get_lambda(self.point2.press.value, self.point1.press.value)
        # lambda_st = self.__get_lambda(self.std.point2.press.value, self.std.point1.press.value)
        # Q0_st = self.Q0 * (qv_st * lambda_st) / qv * lambda_

    def get_points_params(self):
        # Точка 1
        o = self.t['o'] + 273.15
        o = round(o)
        vs = self.t['vs'] + 273.15
        vs = round(vs)
        cond = self.t['cond'] + 273.15
        cond = round(cond)
        p_1 = self.__query('amm_wet', 'P', t=o)[0]  # type: float
        p_2 = self.__query('amm_wet', 'P', t=cond)[0]  # type: float
        p_1_tab_min, p_1_tab_max = self.__get_near_p(p_1)
        p_2_tab_min, p_2_tab_max = self.__get_near_p(p_2)
        v_1, h_1, s_1, p_tab_min, p_tab_max = self.__get_over_p_t(p_1, vs, p_1_tab_min, p_1_tab_max)
        self.point1.temp.value = self.t['vs']
        self.point1.press.value = p_1 / 10
        self.point1.ent.value = h_1
        self.point1.cap.value = v_1
        self.point1.state = 'over'
        self.point1.human_name = '1'
        self.point1.print_all()

        # точка 1'
        v_1d, h_1d = self.__query('amm_wet', 'vdd, hdd', p=p_1)
        self.point1d.temp.value = self.t['o']
        self.point1d.press.value = p_1 / 10
        self.point1d.ent.value = h_1d
        self.point1d.cap.value = v_1d
        self.point1d.state = 'wet'
        self.point1d.human_name = '1\''
        self.point1d.print_all()

        t_2, h_2, v_2 = self.__get_over_p_s(p_2, s_1, p_2_tab_min, p_2_tab_max)

        self.point2.temp.value = t_2 - 273.15
        self.point2.press.value = p_2 / 10
        self.point2.ent.value = h_2
        self.point2.cap.value = v_2
        self.point2.state = 'over'
        self.point2.human_name = '2'
        self.point2.print_all()

        v_2d, h_2d = self.__get_wet_p_max(p_2)
        self.point2d.temp.value = self.t['cond']
        self.point2d.press.value = self.point2.press.value
        self.point2d.ent.value = h_2d
        self.point2d.cap.value = v_2d
        self.point2d.state = 'wet'
        self.point2d.human_name = '2\''
        self.point2d.print_all()

        v_3d, h_3d = self.__get_wet_p_min(p_2)
        self.point3d.temp.value = self.t['cond']
        self.point3d.press.value = self.point2.press.value
        self.point3d.ent.value = h_3d
        self.point3d.cap.value = v_3d
        self.point3d.state = 'wet'
        self.point3d.human_name = '3\''
        self.point3d.print_all()

        h_3, v_3 = self.__get_point3(self.point2d.press.value * 10, self.t['over'] + 273.15)
        self.point3.temp.value = self.t['over']
        self.point3.press.value = self.point2.press.value
        self.point3.ent.value = h_3
        self.point3.cap.value = v_3
        self.point3.state = 'liquid'
        self.point3.human_name = '3'
        self.point3.print_all()

        v_4, v_4d = self.__get_points_4_4d(o, h_3, h_3d)

        self.point4.temp.value = self.t['o']
        self.point4.press.value = self.point1.press.value
        self.point4.ent.value = self.point3.ent.value
        self.point4.cap.value = v_4
        self.point4.state = 'wet'
        self.point4.human_name = '4'
        self.point4.print_all()

        self.point4d.temp.value = self.t['o']
        self.point4d.press.value = self.point1.press.value
        self.point4d.ent.value = self.point3d.ent.value
        self.point4d.cap.value = v_4d
        self.point4d.state = 'wet'
        self.point4d.human_name = '4\''
        self.point4d.print_all()

    def get_main_temps(self):
        # найдем t0 - температуру кипения хладагента в испарителе
        self.t['in_r'] = self.t['out_r'] + 2
        med_t_r = (self.t['in_r'] + self.t['out_r']) / 2
        self.t['o'] = med_t_r - 5

        # найдем температуру входящей воды конденсатора t_in_k
        self.t['in_k'] = self.__get_t_in_k(self.city)
        self.t['out_k'] = self.t['in_k'] + 4
        med_t_k = (self.t['in_k'] + self.t['out_k']) / 2
        self.t['cond'] = med_t_k + 5
        self.t['over'] = self.t['cond'] - 4
        self.t['vs'] = self.t['o'] + 8
        for key in self.t:
            print('{0} is {1}'.format(key, self.t[key]))


main = ColdSystem(115, 'R717', 'Астрахань', -15)
main.get_main_temps()
main.get_points_params()
main.plot_graph()
main.get_std_params()
main.get_charact_params()
