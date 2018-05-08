import sqlite3
import pickle
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
        print('Для данной точки')
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

    def make_cycle(self):
        # найдем параметры 1' точки
        self.point1.temp.value = self.t['vs']
        self.point1.press.value = self.get_param('давление', self.point1.temp.amount, self.point1.temp.value)
        self.point1.ent.value = self.get_param('энтальпия', self.point1.temp.amount, self.point1.temp.value)
        self.point1.cap.value = self.get_param('удельный объем', self.point1.temp.amount, self.point1.temp.value)
        self.point1.state = self.point1.get_state(self.point1)
        self.point1.human_name = '1'
        self.point1.print_all()

        self.point1d.temp.value = self.t['o']
        self.point1d.press.value = self.point1.press.value
        self.point1d.ent.value = self.get_param('энтальпия', self.point1d.temp.amount, self.point1d.temp.value)
        self.point1d.cap.value = self.get_param('удельный объем', self.point1d.temp.amount, self.point1d.temp.value)
        self.point1d.state = self.point1d.get_state(self.point1d)
        self.point1d.human_name = '1\''
        self.point1d.print_all()

        self.point2.temp.value = self.get_param('температура', self.point2.temp.amount, self.point2.temp.value)
        self.point2.press.value = self.get_param('давление', self.point2.temp.amount, self.point2.temp.value)
        self.point2.ent.value = self.get_param('энтальпия', self.point2.temp.amount, self.point2.temp.value)
        self.point2.cap.value = self.get_param('удельный объем', self.point2.temp.amount, self.point2.temp.value)
        self.point2.state = self.point2.get_state(self.point2)
        self.point2.human_name = '2'
        self.point2.print_all()

        self.point2d.temp.value = self.t['over']
        self.point2d.press.value = self.point2.press.value
        self.point2d.ent.value = self.get_param('энтальпия', self.point2d.temp.amount, self.point2d.temp.value)
        self.point2d.cap.value = self.get_param('удельный объем', self.point2d.temp.amount, self.point2d.temp.value)
        self.point2d.state = self.point2d.get_state(self.point2d)
        self.point2d.human_name = '2\''
        self.point2d.print_all()

        self.point3.temp.value = self.point2d.temp.value
        self.point3.press.value = self.point2.press.value
        self.point3.ent.value = self.get_param('энтальпия', self.point3.temp.amount, self.point3.temp.value)
        self.point3.cap.value = self.get_param('удельный объем', self.point3.temp.amount, self.point3.temp.value)
        self.point3.state = self.point3.get_state(self.point3)
        self.point3.human_name = '3'
        self.point3.print_all()

        self.point3d.temp.value = self.get_param('температура', self.point3d.temp.amount, self.point3d.temp.value)
        self.point3d.press.value = self.point2.press.value
        self.point3d.ent.value = self.get_param('энтальпия', self.point3d.temp.amount, self.point3d.temp.value)
        self.point3d.cap.value = self.get_param('удельный объем', self.point3d.temp.amount, self.point3d.temp.value)
        self.point3d.state = self.point3d.get_state(self.point3d)
        self.point3d.human_name = '3\''
        self.point3d.print_all()

        self.point4.temp.value = self.t['o']
        self.point4.press.value = self.point1.press.value
        self.point4.ent.value = self.get_param('энтальпия', self.point4.temp.amount, self.point4.temp.value)
        self.point4.cap.value = self.get_param('удельный объем', self.point4.temp.amount, self.point4.temp.value)
        self.point4.state = self.point4.get_state(self.point4)
        self.point4.human_name = '4'
        self.point4.print_all()

        self.point4d.temp.value = self.t['o']
        self.point4d.press.value = self.point1.press.value
        self.point4d.ent.value = self.get_param('энтальпия', self.point4d.temp.amount, self.point4d.temp.value)
        self.point4d.cap.value = self.get_param('удельный объем', self.point4d.temp.amount, self.point4d.temp.value)
        self.point4d.state = self.point4d.get_state(self.point4d)
        self.point4d.human_name = '4\''
        self.point4d.print_all()

    def make_current_test_cycle(self):
        # найдем параметры 1' точки
        self.point1.temp.value = self.t['vs']
        self.point1.press.value = 0.196
        self.point1.ent.value = 1675
        self.point1.cap.value = 0.61
        self.point1.state = 'over'
        self.point1.human_name = '1'
        self.point1.print_all()

        self.point1d.temp.value = self.t['o']
        self.point1d.press.value = self.point1.press.value
        self.point1d.ent.value = 1660
        self.point1d.cap.value = 0.6
        self.point1d.state = 'wet'
        self.point1d.human_name = '1\''
        self.point1d.print_all()

        self.point2.temp.value = 129
        self.point2.press.value = 1.32
        self.point2.ent.value = 1965
        self.point2.cap.value = 0.14
        self.point2.state = 'over'
        self.point2.human_name = '2'
        self.point2.print_all()

        self.point2d.temp.value = self.t['over']
        self.point2d.press.value = self.point2.press.value
        self.point2d.ent.value = 1710
        self.point2d.cap.value = 0.09
        self.point2d.state = 'wet'
        self.point2d.human_name = '2\''
        self.point2d.print_all()

        self.point3.temp.value = self.point2d.temp.value
        self.point3.press.value = self.point2.press.value
        self.point3.ent.value = 565
        self.point3.cap.value = 0
        self.point3.state = 'liquid'
        self.point3.human_name = '3'
        self.point3.print_all()

        self.point3d.temp.value = 38
        self.point3d.press.value = self.point2.press.value
        self.point3d.ent.value = 590
        self.point3d.cap.value = 0
        self.point3d.state = 'liquid'
        self.point3d.human_name = '3\''
        self.point3d.print_all()

        self.point4.temp.value = self.t['o']
        self.point4.press.value = self.point1.press.value
        self.point4.ent.value = 565
        self.point4.cap.value = 0.105
        self.point4.state = 'wet'
        self.point4.human_name = '4'
        self.point4.print_all()

        self.point4d.temp.value = self.t['o']
        self.point4d.press.value = self.point1.press.value
        self.point4d.ent.value = 590
        self.point4d.cap.value = 0.12
        self.point4d.state = 'wet'
        self.point4d.human_name = '4\''
        self.point4d.print_all()

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
                plt.annotate(s='точка ' + point.human_name, xy=(point.ent.value, point.press.value,),
                             xytext=(point.ent.value, point.press.value+0.03,))
            else:
                plt.annotate(s='точка ' + point.human_name, xy=(point.ent.value, point.press.value,),
                             xytext=(point.ent.value, point.press.value-0.03,))
        plt.plot(x, z)
        plt.xlim(500, 2000)
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


main = ColdSystem(115, 'R717', 'Астрахань', -15)
main.get_main_temps()
main.make_current_test_cycle()
main.plot_graph()
