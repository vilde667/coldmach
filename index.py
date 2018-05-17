import sqlite3
import pickle
import matplotlib.pyplot as plt
import numpy as np
import math


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

    @staticmethod
    def get_t_in_k(city):
        conn = sqlite3.connect('main.db')
        c = conn.cursor()
        t = (city,)
        c.execute('SELECT * FROM citys WHERE city=?', t)
        response = c.fetchone()
        return float(response[1].replace(',', '.'))

    def get_main_temps(self):
        # найдем t0 - температуру кипения хладагента в испарителе
        self.t['in_r'] = self.t['out_r'] + 2
        med_t_r = (self.t['in_r'] + self.t['out_r']) / 2
        self.t['o'] = med_t_r - 5

        # найдем температуру входящей воды конденсатора t_in_k
        self.t['in_k'] = self.get_t_in_k(self.city)
        self.t['out_k'] = self.t['in_k'] + 4
        med_t_k = (self.t['in_k'] + self.t['out_k']) / 2
        self.t['cond'] = med_t_k + 5
        self.t['over'] = self.t['cond'] - 4
        self.t['vs'] = self.t['o'] + 8
        for key in self.t:
            print('{0} is {1}'.format(key, self.t[key]))

    @staticmethod
    def get_param(param, ret_par, ret):
        r = input('введите значение {0} при {1} = {2}\n'.format(param, ret_par, ret))
        return float(r)

    def get_points_params(self):
        #Точка 1
        o = self.t['o'] + 273.15
        o = round(o)
        vs = self.t['vs'] + 273.15
        vs = round(vs)
        cond = self.t['cond'] + 273.15
        cond = round(cond)
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT * FROM amm_wet WHERE T={0};'.format(o))
        req = cur.fetchone()
        p_1 = req[2]
        cur.execute('SELECT * FROM amm_wet WHERE T={0};'.format(cond))
        req = cur.fetchone()
        p_2 = req[2]
        p_1_tab_min, p_1_tab_max = self.get_near_p(p_1)
        p_2_tab_min, p_2_tab_max = self.get_near_p(p_2)
        v_1, h_1, s_1, p_tab_min, p_tab_max = self.get_over_p_t(p_1, vs, p_1_tab_min, p_1_tab_max)
        self.point1.temp.value = self.t['vs']
        self.point1.press.value = p_1 / 10
        self.point1.ent.value = h_1
        self.point1.cap.value = v_1
        self.point1.state = 'over'
        self.point1.human_name = '1'
        self.point1.print_all()

        #точка 1'
        cur.execute('SELECT vdd, hdd FROM amm_wet WHERE P={0};'.format(p_1))
        res_1d = cur.fetchone()
        v_1d, h_1d = res_1d
        self.point1d.temp.value = self.t['o']
        self.point1d.press.value = p_1 / 10
        self.point1d.ent.value = h_1d
        self.point1d.cap.value = v_1d
        self.point1d.state = 'wet'
        self.point1d.human_name = '1\''
        self.point1d.print_all()
        conn.close()

        t_2, h_2, v_2 = self.get_over_p_s(p_2, s_1, p_2_tab_min, p_2_tab_max)

        self.point2.temp.value = t_2 - 273.15
        self.point2.press.value = p_2 / 10
        self.point2.ent.value = h_2
        self.point2.cap.value = v_2
        self.point2.state = 'over'
        self.point2.human_name = '2'
        self.point2.print_all()

        v_2d, h_2d = self.get_wet_p_max(p_2)
        self.point2d.temp.value = self.t['cond']
        self.point2d.press.value = self.point2.press.value
        self.point2d.ent.value = h_2d
        self.point2d.cap.value = v_2d
        self.point2d.state = 'wet'
        self.point2d.human_name = '2\''
        self.point2d.print_all()

        v_3d, h_3d = self.get_wet_p_min(p_2)
        self.point3d.temp.value = self.t['cond']
        self.point3d.press.value = self.point2.press.value
        self.point3d.ent.value = h_3d
        self.point3d.cap.value = v_3d
        self.point3d.state = 'wet'
        self.point3d.human_name = '3\''
        self.point3d.print_all()

        h_3, v_3 = self.get_point3()
        self.point3.temp.value = self.t['over']
        self.point3.press.value = self.point2.press.value
        self.point3.ent.value = h_3
        self.point3.cap.value = v_3
        self.point3.state = 'liquid'
        self.point3.human_name = '3'
        self.point3.print_all()

        v_4, v_4d = self.get_points_4_4d(o, h_3, h_3d)

        self.point4.temp.value = self.t['o']
        self.point4.press.value = self.point1.press.value
        self.point4.ent.value = self.point3.ent.value
        self.point4.cap.value = v_4
        self.point4.state = 'wet'
        self.point4.human_name = '4'
        self.point4.print_all()

        self.point4d.temp.value = self.t['o']
        self.point4d.press.value = self.point1.press.value
        self.point4d.ent.value =  self.point3d.ent.value
        self.point4d.cap.value = v_4d
        self.point4d.state = 'wet'
        self.point4d.human_name = '4\''
        self.point4d.print_all()

    def get_points_4_4d(self, t, h_4, h_4d):
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT hd, hdd, vd, vdd FROM amm_wet WHERE T={0}'.format(t))
        hd, hdd, vdd, vd = cur.fetchone()
        x_4 = (h_4 - hd) / (hdd - hd)
        v_4 = vd + x_4 * (vdd - vd)
        x_4d = (h_4d - hd) / (hdd - hd)
        v_4d = vd + x_4d * (vdd - vd)
        return v_4, v_4d

    def get_point3(self):
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        p = self.point2d.press.value * 10
        t = self.t['over'] + 273.15
        t_tab_min, t_tab_max = self.get_near_t(t)
        p_tab_min, p_tab_max = self.get_near_p(p)
        res_all = []
        cur.execute('SELECT v,h FROM amm_over WHERE P={0} AND T={1}'.format(p_tab_min, int(t_tab_min)))
        res_all.append(cur.fetchone())
        cur.execute('SELECT v,h FROM amm_over WHERE P={0} AND T={1}'.format(p_tab_max, int(t_tab_min)))
        res_all.append(cur.fetchone())
        zipped = tuple(zip(res_all[0], res_all[1]))
        v = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[0][1] - zipped[0][0]) + zipped[0][0]
        h = (t - t_tab_min) / (t_tab_max - t_tab_min) * (zipped[1][1] - zipped[1][0]) + zipped[1][0]
        conn.close()
        return h, v

    def get_wet_p_max(self, p=None):
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT vdd, hdd FROM amm_wet WHERE P={0}'.format(p))
        v, h = cur.fetchone()
        conn.close()
        return v, h

    def get_wet_p_min(self, p=None):
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT vd, hd FROM amm_wet WHERE P={0}'.format(p))
        v, h = cur.fetchone()
        conn.close()
        return v, h

    def get_near_p(self, p):
        delta_p = 10
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT * FROM amm_over WHERE P>{0} AND P<{1} AND T=300'.format(p - delta_p, p + delta_p))
        response = cur.fetchall()
        near_minus = []
        near_plus = []
        for x in response:
            if p - x[1] >= 0:
                near_plus.append(p - x[1])
            else:
                near_minus.append(p - x[1])
        p_tab_min_ = max(near_minus)
        p_tab_max_ = min(near_plus)
        p_tab_min = p - p_tab_min_
        p_tab_max = p - p_tab_max_
        conn.close()
        return p_tab_min, p_tab_max

    def get_near_t(self, t):
        delta_t = 10
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT T FROM amm_over WHERE T>{0} AND T<{1} AND P=3.0'.format(t - delta_t, t + delta_t))
        response = cur.fetchall()
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
        conn.close()
        return t_tab_max, t_tab_min

    def get_over_p_t(self, p=None, t=None, p_tab_min=None, p_tab_max=None):
        print(t, p)
        delta_t = 5  # type: int
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute('SELECT * FROM amm_over WHERE P={0} AND T>{1} AND T<{2}'.format(
            p_tab_min, t - delta_t, t + delta_t))
        response = cur.fetchall()

        t_tab = response[0][0]
        if t_tab - t >= 0:
            t_tab_max = t_tab
            t_tab_min = t_tab - 10
        else:
            t_tab_min = t_tab
            t_tab_max = t_tab + 10

        res_all = []
        cur.execute('SELECT v,h,s FROM amm_over WHERE P={0} AND T={1}'.format(p_tab_min, int(t_tab_min)))
        res_all.append(cur.fetchone())
        cur.execute('SELECT v,h,s FROM amm_over WHERE P={0} AND T={1}'.format(p_tab_max, int(t_tab_min)))
        res_all.append(cur.fetchone())
        cur.execute('SELECT v,h,s FROM amm_over WHERE P={0} AND T={1}'.format(p_tab_min, int(t_tab_max)))
        res_all.append(cur.fetchone())
        cur.execute('SELECT v,h,s FROM amm_over WHERE P={0} AND T={1}'.format(p_tab_max, int(t_tab_max)))
        res_all.append(cur.fetchone())
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
        conn.close()
        return v, h, s, p_tab_min, p_tab_max

    def get_over_p_s(self, p=None, s=None, p_tab_min=None, p_tab_max=None):
        print(s, p)
        delta_s = 0.2
        conn = sqlite3.connect('main.db')
        cur = conn.cursor()
        cur.execute(
            'SELECT t,s FROM amm_over WHERE P={0} AND S>{1} AND S<{2};'.format(p_tab_min, s - delta_s, s + delta_s))
        response = cur.fetchall()
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
        cur.execute(
            'SELECT t,s FROM amm_over WHERE P={0} AND S>{1} AND S<{2}'.format(p_tab_max, s - delta_s, s + delta_s))
        response = cur.fetchall()
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
        res_all = []
        cur.execute(
            'SELECT t,h,v FROM amm_over WHERE P={0} AND S={1};'.format(p_tab_min, s_tab_min1))
        res_all.append(cur.fetchone())
        cur.execute(
            'SELECT t,h,v FROM amm_over WHERE P={0} AND S={1};'.format(p_tab_min, s_tab_max1))
        res_all.append(cur.fetchone())
        cur.execute(
            'SELECT t,h,v FROM amm_over WHERE P={0} AND S={1};'.format(p_tab_max, s_tab_min2))
        res_all.append(cur.fetchone())
        cur.execute(
            'SELECT t,h,v FROM amm_over WHERE P={0} AND S={1};'.format(p_tab_max, s_tab_max2))
        res_all.append(cur.fetchone())
        zipped = tuple(zip(res_all[0], res_all[1], res_all[2] ,res_all[3]))
        t1 = (s - s_tab_min1) / (s_tab_max1 - s_tab_min1) * (zipped[0][1] - zipped[0][0]) + zipped[0][0]
        t2 = (s - s_tab_min2) / (s_tab_max2 - s_tab_min2) * (zipped[0][3] - zipped[0][2]) + zipped[0][2]
        t = (p - p_tab_min) / (p_tab_max - p_tab_min) * (t2 - t1) + t1
        h1 = (s - s_tab_min1) / (s_tab_max1 - s_tab_min1) * (zipped[1][1] - zipped[1][0]) + zipped[1][0]
        h2 = (s - s_tab_min2) / (s_tab_max2 - s_tab_min2) * (zipped[1][3] - zipped[1][2]) + zipped[1][2]
        h = (p - p_tab_min) / (p_tab_max - p_tab_min) * (h2 - h1) + h1
        v1 = (s - s_tab_min1) / (s_tab_max1 - s_tab_min1) * (zipped[2][1] - zipped[2][0]) + zipped[2][0]
        v2 = (s - s_tab_min2) / (s_tab_max2 - s_tab_min2) * (zipped[2][3] - zipped[2][2]) + zipped[2][2]
        v = (p - p_tab_min) / (p_tab_max - p_tab_min) * (v2 - v1) + v1
        conn.close()
        return t, h, v

    def plot_graph(self):
        x = np.arange(self.point1.ent.value, self.point2.ent.value)
        p = self.point2.press.value
        r = self.point2.ent.value
        g = (p - self.point1.press.value)/((self.point1.ent.value - r) ** 2)
        z = -g * (x - r) ** 2 + p
        self.points = [self.point1, self.point1d, self.point2, self.point2d, self.point3, self.point3d, self.point4, self.point4d]
        plt.plot([self.point2.ent.value, self.point3.ent.value], [self.point2.press.value, self.point3.press.value])
        plt.plot([self.point3.ent.value, self.point4.ent.value], [self.point3.press.value, self.point4.press.value])
        plt.plot([self.point4.ent.value, self.point1.ent.value], [self.point4.press.value, self.point1.press.value])
        for point in self.points:
            plt.scatter(point.ent.value, point.press.value)
            if '\'' in point.human_name:
                plt.annotate(s=point.human_name, xy=(point.ent.value, point.press.value,),
                             xytext=(point.ent.value, point.press.value+0.03,))
            else:
                plt.annotate(s=point.human_name, xy=(point.ent.value, point.press.value,),
                             xytext=(point.ent.value + 10, point.press.value-0.06,))
        plt.plot(x, z)
        plt.xlim(200, 2100)
        plt.ylim(0.1, 1.5)
        plt.show()

    def try_to_pickle(self):
        pickler = [self.point1, self.point1d, self.point2, self.point2d, self.point3, self.point3d, self.point4,
                   self.point4d]
        with open('pickler.pkl', 'wb') as file:
            pickle.dump(pickler, file)
            print('pickling good!')

    def try_to_unpickle(self):
        with open('pickler.pkl', 'rb') as file:
            pickler = pickle.load(file)
            self.point1 = pickler[0]
            self.point1d = pickler[1]
            self.point2 = pickler[2]
            self.point2d = pickler[3]
            self.point3 = pickler[4]
            self.point3d = pickler[5]
            self.point4 = pickler[6]
            self.point4d = pickler[7]


main = ColdSystem(115, 'R717', 'Салехард', -15)
main.get_main_temps()
main.get_points_params()
main.plot_graph()
