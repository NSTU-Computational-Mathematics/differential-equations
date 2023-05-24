import math

DIMENSION = 4

class AdamsMethod:
    points = []

    def __init__(self, first_points, last_abscissa, factors, h, func):
        self.points = first_points
        self.func = func
        self.factors = factors
        self.last_values = [func(point[0], point[1]) for point in first_points]
        self.h = h
        self.last_abscissa = last_abscissa

    def calc_points(self):
        while True:
            point = self.__next()
            self.points.append(point)

            if point[0] >= self.last_abscissa:
                break;

        dots_x = []
        dots_y = []

        for point in self.points:
            dots_x.append(point[0])
            dots_y.append(point[1])
        
        return dots_x, dots_y

    def __next(self):
        last_ordinate = self.points[-1][1];
        last_abscissa = self.points[-1][0];

        last_values_num = len(self.last_values)

        next_abscissa = last_abscissa + self.h
        next_ordinate = last_ordinate + (self.h / 24) * \
            sum([self.last_values[i] * self.factors[i] for i in range(last_values_num)])

        del self.last_values[0]
        self.last_values.append(self.func(next_abscissa, next_ordinate))

        return (next_abscissa, next_ordinate)
