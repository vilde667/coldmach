import sqlite3

conn = sqlite3.connect('main.db')
c = conn.cursor()

p = []

#c.execute('SELECT P,T,v,h,s,cp,k,a FROM amm_water ORDER BY P ASC')
c.execute('SELECT P FROM amm_over ORDER BY P ASC')

res = c.fetchall()
for x in res:
    p.append(x[0])

p = list(set(p))

p.sort()

for el in p:
    c.execute('SELECT P,T,v,h,s,cp,k,a FROM amm_over WHERE P={0} ORDER BY T ASC'.format(el))

    res = c.fetchall()


    if len(res) > 1:
        for line in range(0, len(res)-1):
            current = res[line]
            nxt = res[line + 1]

            print(current)
            print(nxt) 

            for n in range(1, 10):
                print(current[0])
                print(current[1] + (nxt[1]-current[1])/10*n)
                print(current[2] + (nxt[2]-current[2])/10*n)
                print(current[3] + (nxt[3]-current[3])/10*n)
                print(current[4] + (nxt[4]-current[4])/10*n)
                print(current[5] + (nxt[5]-current[5])/10*n)
                print(current[6] + (nxt[6]-current[6])/10*n)
                print(current[7] + (nxt[7]-current[7])/10*n)

                tup = (current[0], current[1] + (nxt[1]-current[1])/10*n, current[2] + (nxt[2]-current[2])/10*n, 
                    current[3] + (nxt[3]-current[3])/10*n, current[4] + (nxt[4]-current[4])/10*n,
                    current[5] + (nxt[5]-current[5])/10*n, current[6] + (nxt[6]-current[6])/10*n,
                    current[7] + (nxt[7]-current[7])/10*n)

                c.execute('INSERT INTO amm_over(P,T,v,h,s,cp,k,a) VALUES (?,?,?,?,?,?,?,?)', tup)

                conn.commit()