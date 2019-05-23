import matplotlib.pyplot as plt
import math

sections = ((-2.5, -1.6), (-1.4, -0.84), (1.5, 2.27))

def f(x):
    return ((x ** 9  + math.pi) * math.cos(math.log(x ** 2 + 1))) / math.exp(x ** 2) - x / 2018

def draw_function(f, a, b, step):

    x = list()
    y = list()

    i = a
    while i <= b:
        x.append(i)
        y.append(f(i))

        i += step

    x_ticks = list()
    i = a
    while i <= b:
        x_ticks.append(i)

        i += 1

    print((x[0], y[0]))
    print((x[-1], y[-1]))

    plt.plot(x, y)
    plt.grid(True)
    plt.xticks(x_ticks)

    plt.show()

def sign(x):
    if x < 0:
        return -1
    elif x > 0:
        return 1

    return 0

def bisection(sections, f, section_length):
    result_sections = [list(section) for section in sections]

    for section in result_sections:
        while section[1] - section[0] > section_length:
            mid = (section[1] + section[0]) / 2
            if sign(f(mid)) != sign(f(section[0])):
                section[1] = mid
            else:
                section[0] = mid

    return result_sections

def newton_discrete(sections, f, h, eps):
    pass

def newton(sections, f, ):
    pass

if __name__ == "__main__":
    # a = -27 - math range error

    # draw_function(f, -26, 20, 0.01)

    result_sections = bisection(sections, f, 10e-4)




