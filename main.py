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


euler_method(0, 1, 3, 0.000001)
runge_kutta_merson_method(0, 1, 3, 0.01, 0.01)
