'''Лабораторная работа №3
Моделирование расчёта толщины окисла
от времени (30 мин). Толщина первого
слоя 0,1 мкм (уровень легирования 10^15 см^-3).
n-Si 10,10^17 см^-3.
Режим окисления O2, T = 1373 K, P = 2 атм.
'''

import math
import matplotlib.pyplot as plt

t = 30
T = 1373
P = 2
ni = 1e17
n = 1e15
xi = 0.1 #мкм
k = 8.62e-5

#Посчитаем зависимость линейной константы от концентрации легирующей примеси
Eg = 1.12 - (4.73e-4) * T ** 2 / (T + 636)
Ei = Eg / 2 - k * T / 4
Cm = math.exp((Ei + 0.57 - Eg) / (k * T))
Cmm = math.exp((2 * Ei + 1.25 - 3 * Eg) / (k * T))
Cp = math.exp((0.35 - Ei) / (k * T))
Vn = (1 + Cp * ni / n) + Cm * n / ni + Cmm * ((n / ni) ** 2)/(1 + Cp + Cm + Cmm)
gammal = 2620 * math.exp(-1.1 / k / T)

#Посчитаем зависимость параболической константы от концентрации легирующей примеси
gammap = (9.63e-16) * math.exp(2.83 / k / T)

A0 = 12.9
EA = 1.23
B = A0 * math.exp(-EA / (k * T)) * (1 + gammap * (n ** 0.22))
A0 = 1.04e5
EA = 2
BA = A0 * math.exp(-EA / (k * T)) * (P ** 0.75) * (1 + gammal * (Vn - 1))

x = []
#Найдём решения квадратного уравнения
for i in range(0, t + 1):
    a = 1/B
    b = 1/BA
    c = -(xi ** 2)/B - xi/BA - i
    dis = b ** 2 - 4 * a * c
    x.append((-b + math.sqrt(dis)) / (2 * a))

plt.suptitle('Зависимость толщины окисла от времени.\n Финальная толщина = {0:.3f} мкм'.format(x[t]))
plt.axis([0, t, xi, x[t]])
t = [i for i in range(0,t + 1)]
plt.plot(t, x)
plt.xlabel('Время, мин')
plt.ylabel('Глубина, мкм')
plt.show()
