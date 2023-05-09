from peewee import *

db = SqliteDatabase('./db/main.db')


def get_approx_value(table, param, param_str, known_param, known_param_str, known_value):
    known_values = []
    results = table.select()
    for r in results:
        known_values.append(getattr(r, known_param_str))
    res = list(set(known_values))
    res.sort()

    known_value_min = min(res, key=lambda x: abs(x - known_value))

    if known_value_min > known_value:
        known_value_min = res[res.index(known_value_min) - 1]

    known_value_max = res[res.index(known_value_min) + 1]

    if known_value == known_value_min:
        results = table.select().where(known_param == known_value_min)
        for r in results:
            return getattr(r, param_str)
    elif known_value == known_value_max:
        results = table.select().where(known_param == known_value_max)
        for r in results:
            return getattr(r, param_str)
    else:
        param_min = 0
        param_max = 0
        results = table.select().where(known_param == known_value_min)

        for r in results:
            param_min = getattr(r, param_str)

        results = table.select().where(known_param == known_value_max)

        for r in results:
            param_max = getattr(r, param_str)

        k = (known_value - known_value_min) / (known_value_max - known_value_min)
        return param_min + (param_max - param_min) * k


class Cities(Model):
    __doc__ = '''Список городов со средими температурами воздеха в летний период'''

    city = CharField()
    t_air_out = IntegerField()

    class Meta:
        database = db


class FeedRatio(Model):
    __doc__ = '''Коэффициент подачи для реального компреммора для разных хладагентов'''

    compression_ratio = IntegerField()
    R717 = FloatField()
    R22 = FloatField()
    R12 = FloatField()

    class Meta:
        database = db


class Compressors(Model):
    __doc__ = '''Список компрессоров с физическими параметрами и значениями холодопроизводительности 
    и потребляемой мощности на разных хладагентах'''
    mark = CharField()
    diameter = FloatField()
    piston_stroke = FloatField()
    theor_v = FloatField()
    frequency = IntegerField()
    Q_R22 = FloatField()
    N_R22 = FloatField()
    Q_R12 = FloatField()
    N_R12 = FloatField()
    Q_R717 = FloatField(null=True)
    N_R717 = FloatField(null=True)

    class Meta:
        database = db


class R717Liquid(Model):
    __doc__ = """Значения теплофизических и термодинамических параметров фреона R717(аммиак) в области жидкости"""

    T = FloatField()
    P = FloatField()
    v = FloatField()
    h = FloatField()
    s = FloatField()
    cp = FloatField()
    k = FloatField()
    a = FloatField()

    class Meta:
        database = db


class R717Wet(Model):
    __doc__ = """Значения теплофизических и термодинамических параметров фреона R717(аммиак) 
    в области влажного насыщенного пара"""

    T = FloatField()
    P = FloatField()
    vd = FloatField()
    vdd = FloatField()
    rod = FloatField()
    rodd = FloatField()
    hd = FloatField()
    hdd = FloatField()
    r = FloatField()
    sd = FloatField()
    sdd = FloatField()

    class Meta:
        database = db


class R717Over(Model):
    __doc__ = """Значения теплофизических и термодинамических параметров 
    фреона R717(аммиак) в области перегретого пара"""

    T = FloatField()
    P = FloatField()
    v = FloatField()
    h = FloatField()
    s = FloatField()
    cp = FloatField()
    k = FloatField()
    a = FloatField()

    class Meta:
        database = db
