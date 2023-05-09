import sqlite3

conn = sqlite3.connect('main.db')
c = conn.cursor()

c.execute('INSERT INTO amm_over(P,T,v,h,s,cp,k,a) VALUES (?,?,?,?,?,?,?,?)', (280.0, 370.0, 0.01932, 865.0, 5.3, 4.94, 21.0, 1080.0))

conn.commit()