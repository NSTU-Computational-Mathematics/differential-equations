def f(x, y):
    return y - x ** 2


def euler_method(x0, y0, xn, n):
    h = (xn - x0) / n

    xn = x0
    yn = y0

    for i in range(n):
        yn1 = yn + h * f(xn, yn)
        xn1 = xn + h
        xn, yn = xn1, yn1

    print("x =", xn, "y =", yn)


def runge_kutta_merson_method(x0, y0, xn, n):
    h = (xn - x0) / n

    xn = x0
    yn = y0

    for i in range(n):
        k1 = f(xn, yn)
        k2 = f(xn + h / 2, yn + h / 2 * k1)
        k3 = f(xn + h / 2, yn + h / 2 * k2)
        k4 = f(xn + h, yn + h * k3)
        yn1 = yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        xn1 = xn + h
        xn, yn = xn1, yn1

    print("x =", xn, "y =", yn)


euler_method(0, 1, 3, 400)
runge_kutta_merson_method(0, 1, 3, 30)
