from peewee import *

db = SqliteDatabase('./db/main.db')


class Cities(Model):
    __doc__ = '''Список городов со средими температурами воздеха в летний период'''

    city = CharField()
    t_air_out = IntegerField()

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
