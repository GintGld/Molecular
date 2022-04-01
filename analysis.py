import matplotlib.pyplot as plt
import numpy as np
import modin.pandas as md

bin = 0.03

def mnk(a, mi, ma):
    xy = 0
    xx = 0
    for key in a:
        #if (float(key) >= 0.5 * mi and float(key) <= 0.5 * ma):
        xy += np.log(a[key] / a[0]) * key ** 2
        xx += np.log(a[key] / a[0]) ** 2
    k = xy / xx
    return k

def energy():
    en = open('Energy.txt', 'r')
    s = en.readline()
    s = s.rstrip()
    dt =float(s)
    s = en.readline()
    s = s.rstrip()
    k = []
    p = []
    f = []
    while(s != 'end'):
        s = s.split(',')
        k.append(float(s[0]))
        p.append(float(s[1]))
        f.append(float(s[0]) + float(s[1]))
        s = en.readline()
        s = s.rstrip()
    en.close()
    x = np.arange(0, len(k) * dt, dt)
    plt.plot(x, k, label='potential energy')
    plt.plot(x, p, label='kinetic energy')
    plt.plot(x, f, label='full energy')
    plt.legend(loc='best', fontsize=12)
    plt.show()

def velocity():
    vel = open('velocity.txt', 'r')
    s = vel.readline()
    s = s.rstrip()
    x = dict()
    y = dict()
    z = dict()
    v = dict()
    m = dict()
    m['max_x'] = 0
    m['max_y'] = 0
    m['max_z'] = 0
    m['min_x'] = 0
    m['min_y'] = 0
    m['min_z'] = 0
    m['max_v'] = 0
    counter = 0
    while(s != 'end'):
        s = s.split(',')
        x_0, y_0, z_0 = float(s[0]) - (float(s[0]) % bin), float(s[1]) - (float(s[1]) % bin), float(s[2]) - (float(s[2]) % bin)
        r = np.sqrt(float(s[0]) ** 2 + float(s[1]) ** 2 + float(s[2]) ** 2)
        v_0 = r - (r % bin)
        if x_0 in x: x[x_0] += 1
        else: x[x_0] = 1
        if y_0 in y: y[y_0] += 1
        else: y[y_0] = 1
        if z_0 in z: z[z_0] += 1
        else: z[z_0] = 1
        if v_0 in v: v[v_0] += 1
        else: v[v_0] = 1
        m['max_x'] = max(m['max_x'], x_0)
        m['min_x'] = min(m['min_x'], x_0)
        m['max_y'] = max(m['max_y'], y_0)
        m['min_y'] = min(m['min_y'], y_0)
        m['max_z'] = max(m['max_z'], z_0)
        m['min_z'] = min(m['min_z'], z_0)
        m['max_v'] = max(m['max_v'], v_0)
        s = vel.readline()
        s = s.rstrip()
        counter += 1
    vel.close()

    for key in x: x[key] /= counter * bin
    for key in y: y[key] /= counter * bin
    for key in z: z[key] /= counter * bin
    for key in v: v[key] /= counter * bin
    x_list = np.arange(m['min_x'], m['max_x'] + 0.05, 0.05)
    y_list = np.arange(m['min_y'], m['max_y'] + 0.05, 0.05)
    z_list = np.arange(m['min_z'], m['max_z'] + 0.05, 0.05)
    v_list = np.arange(0, m['max_v'] + 0.05, 0.05)
    k = [0] * 3
    k[0] = mnk(x, m['min_x'], m['max_x'])
    k[1] = mnk(y, m['min_y'], m['max_y'])
    k[2] = mnk(z, m['min_z'], m['max_z'])
    k[0] = 1 / k[0]
    k[1] = 1 / k[1]
    k[2] = 1 / k[2]
    t = np.sqrt((k[0] ** 2 + k[1] ** 2 + k[2] ** 2) / 3)

    plt.figure(figsize=(12, 7))
    sp = plt.subplot(2,3,1)
    plt.title('$X$, распределение', fontsize=14)
    plt.plot(x_list, np.sqrt(t / np.pi) * np.exp(-t * x_list ** 2))
    for key in x:
        plt.bar(float(key), float(x[key]), color='blue', width=bin)
        #plt.scatter(float(key), float(x[key]), c='blue')
    
    sp = plt.subplot(2,3,2)
    plt.title('$Y$, распределение', fontsize=14)
    plt.plot(y_list, np.sqrt(t / np.pi) * np.exp(-t * y_list ** 2))
    for key in y:
        plt.bar(float(key), float(y[key]), color='blue', width=bin)
        #plt.scatter(float(key), float(y[key]), c='blue')
    sp = plt.subplot(2,3,3)
    plt.title('$Z$, распределение', fontsize=14)
    plt.plot(z_list, np.sqrt(t / np.pi) * np.exp(-t * z_list ** 2))
    for key in z:
        plt.bar(float(key), float(z[key]), color='blue', width=bin)
        #plt.scatter(float(key), float(z[key]), c='blue')
    sp = plt.subplot(2,3,4)
    plt.title('$X$, линеаризованный', fontsize=14)
    plt.plot(x_list ** 2, -t * x_list ** 2)
    for key in x:
        plt.scatter(float(key) ** 2, np.log(float(x[key]) / x[0]), c='blue', s=3)
    sp = plt.subplot(2,3,5)
    plt.title('$Y$, линеаризованный', fontsize=14)
    plt.plot(y_list ** 2, -t * y_list ** 2)
    for key in y:
        plt.scatter(float(key) ** 2, np.log(float(y[key]) / y[0]), c='blue', s=3)
    sp = plt.subplot(2,3,6)
    plt.title('$Z$, линеаризованный', fontsize=14)
    plt.plot(z_list ** 2, -t * z_list ** 2)
    for key in z:
        plt.scatter(float(key) ** 2, np.log(float(z[key]) / z[0]), c='blue', s=3)

    plt.savefig('velocity\\axes.png')

    plt.figure(figsize=(18, 11))
    sp = plt.subplots()
    plt.title('$V$', fontsize=14)
    plt.plot(v_list, 4 * np.pi * v_list ** 2 * (t / np.pi) ** (1.5) * np.exp(-t * v_list ** 2), label='T = ' + str(round(0.5 / t, 2)))
    for key in v:
        plt.bar(float(key), float(v[key]), color='blue', width=bin)
        #plt.scatter(float(key), float(v[key]), c='blue')
    plt.legend(loc='best', fontsize=12)

    plt.savefig('velocity\\velocity.png')
    #plt.show()
    print('done')
    return

def diffusion():
    d = md.read_csv('diffusion.txt', sep=',')
    d.head()
    print(d)
    return
    N = d.N[0]
    dt = d.dt[0]
    t = dict()
    for step in range(1, 2, 1):
        t[step * dt] = 0
        counter = 0
        #print(t[step * dt])
        for i in range(step, len(d), 1):
            f = 0
            for k in range(N):
                f += (d[str(k) + 'x'][i] - d[str(k) + 'x'][i - step]) ** 2 
                f += (d[str(k) + 'y'][i] - d[str(k) + 'y'][i - step]) ** 2 
                f += (d[str(k) + 'z'][i] - d[str(k) + 'z'][i - step]) ** 2
            t[step * dt] += f / N
            counter += 1
        t[step * dt] /= counter
        print(t[step * dt])
    #plt.scatter(t.keys(), t.values(), c='blue', s=3)
    #plt.show()
    '''
    diff = open('diffusion.txt', 'r')
    s = diff.readline()
    s = s.split(',')
    n = 3 * int(s[0])
    dt = float(s[1])
    N = int(s[2])
    diff.close()
    d = dict()
    for step in range(2, 3, 1): # Рассматриваю различные шаги усреднения
        diff = open('diffusion.txt', 'r')
        s = diff.readline() # Костыль, который обходит первую строку с константами
        d[step * dt] = 0
        prev = [[] * n] * step # Массив, в котором будет храниться предыдущие step координат частиц
        now = [[] * n] * step # Массив, в котором будет храниться Нынешние step координат частиц
        for i in range(0, step, 1): # Считываю начальные условия в prev
            s = diff.readline()
            s = s.rstrip()
            prev[i] = s.split(',')
        for i in range(1, int(N / step), 1): # обход по следующим step строкам
            for g in range(0, step, 1): # Считываю актуальное значение в now
                s = diff.readline()
                s = s.rstrip()
                now[g] = s.split(',')
            f = 0
            for k in range(n): f += float(now[i][k]) ** 2 - float(prev[i][k]) ** 2
            f /= n
            d[step * dt] += f
            prev = now
        d[step * dt] /= int(N / step)
        diff.close()
        print(d[step * dt])
    '''

if __name__ == "__main__":
    #energy()
    #velocity()
    diffusion()