import sqlite3

conn = sqlite3.connect('main.db')
c = conn.cursor()

c.execute('SELECT P,T,v,h,s,cp,k,a FROM amm_over ORDER BY P ASC')

res = c.fetchall()
for x in res:
    for el in x:
        if type(el) == str:
            print(x)