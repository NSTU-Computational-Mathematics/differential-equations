import numpy as np
import matplotlib.pyplot as plt


def f(x, y):
    return y - x ** 2


def euler_method(x0, y0, xn, h):
    n = int((xn - x0) / h)
    xn = x0
    yn = y0
    for i in range(n):
        yn1 = yn + h * f(xn, yn)
        xn1 = xn + h
        xn, yn = xn1, yn1
    print("x =", xn, "y =", yn)
    return xn, yn


def runge_kutta_merson_method(x0, y0, xn, h, eps):
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
    print("x =", x0, "y =", y1)
    return x0, y1


def runge_kutta_4_method(x0, y0, xn, h, eps):
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
    print("x =", x0, "y =", y1)
    return x0, y1


def plotting(dot_x, dot_y, title_name):
    x0, y0 = 0, 1
    xn = 3
    h = 0.01

    x_arr = np.arange(x0, xn + h, h)
    y_arr = np.zeros_like(x_arr)
    y_arr[0] = y0

    for i in range(len(x_arr) - 1):
        k1 = f(x_arr[i], y_arr[i])
        k2 = f(x_arr[i] + h / 2, y_arr[i] + h * k1 / 2)
        k3 = f(x_arr[i] + h / 2, y_arr[i] + h * k2 / 2)
        k4 = f(x_arr[i] + h, y_arr[i] + h * k3)
        y_arr[i + 1] = y_arr[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    plt.plot(x_arr, y_arr)
    plt.plot(dot_x, dot_y, 'o')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title_name)
    plt.grid(True)
    plt.show()


x1, y1 = euler_method(0, 1, 3, 0.000001)
x2, y2 = runge_kutta_merson_method(0, 1, 3, 0.01, 0.01)
x3, y3 = runge_kutta_4_method(0, 1, 3, 0.01, 0.01)

plotting(x1, y1, "Решение дифференциального уравнения y' = y - x^2, метод Эйлера")
plotting(x2, y2, "Решение дифференциального уравнения y' = y - x^2, метод Рунге–Кутты–Мерсона")
plotting(x3, y3, "Решение дифференциального уравнения y' = y - x^2, метод Рунге–Кутты 4-го порядка")
