import numpy as np
import matplotlib.pyplot as plt
from adams_method import AdamsMethod 

def f(x, y):
    return y - x ** 2

def euler_method(x0, y0, xn, h):
    dots_x = []
    dots_y = []
    
    n = int((xn - x0) / h)
    xn = x0
    yn = y0
    for i in range(n):
        yn1 = yn + h * f(xn, yn)
        xn1 = xn + h

        dots_x.append(xn1)
        dots_y.append(yn1)

        xn, yn = xn1, yn1

    return dots_x, dots_y


def runge_kutta_merson_method(x0, y0, xn, h, eps):
    dots_x = []
    dots_y = []

    while x0 < xn:
        k1 = f(x0, y0)
        k2 = f(x0 + h / 3, y0 + h * k1 / 3)
        k3 = f(x0 + h / 3, y0 + h * (k1 + k2) / 6)
        k4 = f(x0 + h / 2, y0 + h * (k1 + 3 * k3) / 8)
        k5 = f(x0 + h, y0 + h * (k1 - 3 * k3 + 4 * k4) / 2)
        y1 = y0 + h * (k1 + 4 * k4 + k5) / 6
        err = abs(y1 - y0) / (1 + abs(y1))
        if err < eps:
            x0 += h
            y0 = y1
            if abs(x0 - xn) < 1e-14:
                break
            h *= min(5, max(0.1, 0.9 * (eps / err) ** 0.2))
        else:
            h *= max(0.1, 0.9 * (eps / err) ** 0.25)
            
        dots_x.append(x0)
        dots_y.append(y1)

    return dots_x, dots_y


def runge_kutta_4_method(x0, y0, xn, h, eps):
    dots_x = []
    dots_y = []

    while x0 < xn:
        k1 = f(x0, y0)
        k2 = f(x0 + h/2, y0 + h*k1/2)
        k3 = f(x0 + h/2, y0 + h*k2/2)
        k4 = f(x0 + h, y0 + h*k3)
        y1 = y0 + h/6 * (k1 + 2*k2 + 2*k3 + k4)
        err = abs(y1 - y0) / (1 + abs(y1))
        if err < eps:
            x0 += h
            y0 = y1
            if abs(x0 - xn) < 1e-14:
                break
            h *= min(5, max(0.1, 0.9 * (eps / err) ** 0.2))
        else:
            h *= max(0.1, 0.9 * (eps / err) ** 0.25)

        dots_x.append(x0)
        dots_y.append(y1)

    return dots_x, dots_y

def adams_method():
    dots_x_3, dots_y_3 = runge_kutta_4_method(0, 1, 4, 0.001, 0.001)
    factors = [-9.0,  37.0, -59.0, 55.0]
    first_points = []
    for i in range(4):
        first_points.append((dots_x_3[i], dots_y_3[i]))
    adams_method = AdamsMethod(first_points, 4, factors, 0.001, f)
    dots_x_4, dots_y_4 = adams_method.calc_points()
    return dots_x_4, dots_y_4

def plotting(dots_x, dots_y, title_name):
    plt.plot(dots_x, dots_y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title_name)
    plt.grid(True)
    plt.show()


dots_x, dots_y = adams_method()
plotting(dots_x, dots_y, "")
