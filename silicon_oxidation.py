"""Лабораторная работа №1
1Провести моделирование одномерного распределения
примеси в кремнии послетехнологических процессов
загонки и диффузионной разгонки.
2. Найти положение глубины залегания pn-перехода
 после проведения обоих процессов.
3. Определить режимы проведения диффузионных
 процессов для заданной примеси для
получения заданной глубины pn-перехода.
Подложка - Si<100> p-типа, 2*10^15см^-3, As,
Режим загонки - имплантация. E=150 КэВ, Q=2*10^14 см^-2
Режим разгонки - диффузия. T = 1000 C, t = 15 мин,
Примесь - мышьяк, глубина pn = 0,6 мкм"""

import math
import matplotlib.pyplot as plt


E = 150000 # Энергия имплантации, Эв
Q = 2e14 # Доза имплантации примеси, см^-3
N = 2e15 # количество атомов в мишени, см^-3
p = 0.00005 # Глубина, см
y = []
x = []
T = 1273 # Температура, К
m = 1.67e-27  # Масса атома
e = 1.6e-19  # Элементарный заряд
M1 = 28.086 * 1.67e-27 # Атомная масса мишени
M2 = 74.92  * 1.67e-27# Атомная масса примеси
Z1 = 14  # Заряд ядра мишени
Z2 = 32  # Заряд ядра примеми
a0 = 0.534 # Постоянная решётки мишени, нм
Rp = 0 # Эв * нм
Csi = 0
dRp1 = 0 # Эв * нм
dE = 1
E += 1
t = 900 # Время загонки, сек
C0 = Q
aa = (0.8854 * a0) / ((Z1 ** (2 / 3) + Z2 ** (2 / 3)) ** 0.5) # Постоянная экранирования
pn = []

def Sn(E): # Ядерная тормозная способность
    a = 1.1383
    b = 0.01321
    c = 0.21226
    d = 0.19593
    eps = (aa * M2 * m * E) / (Z1 * Z2 * 14.4 * (M1 + M2))
    if (eps > 10):
        Seps = 0.5 * math.log(eps) / eps
    else:
        Seps = 0.5 * (math.log(1 + a * eps)) / (eps + b * (eps ** c) + d * (eps ** 0.5))
    S = (8.462e-15 * Z1 * Z2 * M1 * Seps * N) / ((M1 + M2) * (Z1 ** 0.23 + Z2 ** 0.23)) # Эв*см^2
    return S


def Se(E): # Электронная тормозная способность
    k = Z1 ** (1 / 6) * 0.793 * (Z1 ** 0.5) * ((M1 + M2) ** 1.5)
    k = k / ((Z2 ** (2 / 3) + Z1 ** (2 / 3)) ** (3 / 4) * (M1 ** 1.5) * (M2 ** 0.5))
    Cr = 4 * math.pi * (aa ** 2) * M1 * M2 / ((M1 + M2) ** 2)
    Ce = 4 * math.pi * aa * M2 / (Z1 * Z2 * 14.4 * (M1 + M2))
    S = k * (Cr) * (E ** 0.5) * N / (Ce ** 0.5)  # Эв * см^2
    return S


for i in range(1, E): # Расчёт проецированного пробега и его разброса
    Rp = Rp * (1 - M2 * Sn(i) * dE / (2 * M1 * (Se(i) + Sn(i)) * i)) + (dE / (Se(i) + Sn(i)))
    Csi = Csi + 2 * Rp * dE / (Se(i) + Sn(i))
    dRp1 = (dRp1 ** 2 + (Csi - 2 * dRp1 ** 2) * M2 * Sn(i) * dE / (M1 * (Sn(i) + Se(i)) * i))
    dRp = ((Csi - Rp) ** 0.5 - dRp1 ** 2) ** 0.5
    dRp = abs(dRp)
Rp = Rp * 10000000 # Эв * см
dRp = dRp * 1000000 # Эв * см


def C(x): # Распределение Гаусса
    C = (Q / (2.5 * dRp)) * math.exp((x - Rp) ** 2 / (-2 * (Rp ** 2)))
    return C


n = 500 # Точки дискретизации
dx = p/n # Шаг по координате
X = 0 # Координата
i = 0
max_imp = 0
for i in range(0, n): # Запись набора точек профиля распредения в список
    x.append(X)
    y.append(C(X))
    if C(X) > max_imp:
        max_imp = C(X)
    X += dx


def Diff(C): # Расчёт коэффициента диффузии
    k = 0.000086
    Nc = 4.831E+15 * ((1.08 * T) ** 1.5)
    Nv = 4.831E+15 * ((0.59 * T) ** 1.5)
    Eg = 1.21 - (0.00028) * T
    ni = ((Nc * Nv / 2) ** 0.5) * math.exp(-Eg / (2 * k * T))
    Diff = 0.66 * math.exp(-3.44 / (k * T)) + 12 * math.exp (N/ni) * math.exp(-4.05 / (k * T))
    return Diff


# Процесс диффузионной разгонки
b = [0] * (n)
d = [0] * (n)
a = [0] * (n)
r = [0] * (n)
delta = [0] * (n)
lam = [0] * (n)
C = [0] * (n)
#Коэффициенты трехточечного уравнения для левой границы:
b[0] = 0
d[0] = 0
a[0] = 1
r[0] = C0
dt = 1
#Прогоночные коэффициенты из условия b1 = 0:
delta[0] = (-1) * d[0] / a[0]
lam[0] = r[0] / a[0]
#Коэффициенты трехточечного уравнения для правой границы:
b[n - 1] = 0
d[n - 1] = 0
a[n - 1] = 1
r[n - 1] = 0
max_dif = 0
for j in range (1, t):
    for i in range(1, n - 1):
        b[i] = 1
        d[i] = 1
        a[i] = -(2 + ((dx ** 2) / (Diff(C[i]) * dt)))
        r[i] = -((dx ** 2) * C[i]) / (Diff(C[i]) * dt)
    for i in range(1, n - 1):
        delta[i] = -d[i] / (a[i] + b[i] * delta[i - 1])
        lam[i] = (r[i] - b[i] * lam[i - 1]) / (a[i] + b[i] * delta[i - 1])
    C[n - 1] = lam[n - 1]
    for i in range(n - 2, 0, -1):
        C[i] = delta[i] * C[i + 1] + lam[i]
        C[i] = (C[i])
        if C[i] > max_dif:
            max_dif = C[i]
    for i in range(0, n):
        pn.append(abs(Q - C[i]))
C[0] = C[1]


# Расчёт глубины залегания p-n перехода
Min = pn[0]
imin = 0
for i in range(1, n):
    if (pn[i] < pn[i-1]):
        Min = pn[i]
        imin = i


plt.suptitle('Распределение примеси. Глубина p-n перехода = {0:.8f} см'.format(x[imin]))
plt.subplot(121)
plt.title("Загонка")
plt.plot(x, y)
plt.xlabel('Глубина, см')
plt.ylabel('Концентрация, см^-3')
plt.axis([0, p, 0, max_imp])
plt.subplot(122)
plt.title('Разгонка')
plt.plot(x, C)
plt.xlabel('Глубина, см')
plt.ylabel('Концентрация, см^-3')
plt.axis([0, p, 0, max_dif])
plt.tight_layout()
plt.show()
