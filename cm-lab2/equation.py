import matplotlib.pyplot as plt
import math

sections = ((-3.5, -1.6), (-1.4, -0.84), (1.5, 2.27))

BISECTION_SECTION_LENGTH = 1e-4
DISCRETE_NEWTON_H = 1e-5
DISCRETE_NEWTON_EPS = 1e-13
NEWTON_EPS = 1e-15

def f(x):
    return ((x ** 9  + math.pi) * math.cos(math.log(x ** 2 + 1))) / math.exp(x ** 2) - x / 2018

def derivative_f(x):
    return x * math.exp(-x ** 2) * (math.cos(math.log(x ** 2 + 1)) * (9 * x ** 7 - 2 * (x ** 9 + math.pi)) - 2 * (x ** 9 + math.pi) * math.sin(math.log(x ** 2 + 1)) * (x ** 2 + 1) ** -1) - 1 / 2018

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

    plt.plot(x, y)
    plt.grid(True)
    plt.xticks(x_ticks)

    plt.show()

def sign(x):
    return (x >= 0) - (x < 0)

def bisection(sections, f, section_length):
    result_sections = [list(section) for section in sections]
    section_counts = [0] * len(result_sections)

    for i, section in enumerate(result_sections):
        while section[1] - section[0] > section_length:
            mid = (section[1] + section[0]) / 2
            if sign(f(mid)) != sign(f(section[0])):
                section[1] = mid
            else:
                section[0] = mid

            section_counts[i] += 1

    return (result_sections, section_counts)

def discrete_newton(sections, f, h, eps):
    roots = list()
    root_counts = [0] * len(sections)

    for i, section in enumerate(sections):
        prev_x = section[1]

        while True:
            root_counts[i] += 1

            x = prev_x - (h * f(prev_x)) / (f(prev_x + h) - f(prev_x))

            if math.fabs(x - prev_x) < eps:
                break

            prev_x = x

        roots.append(x)

    return (roots, root_counts)

def newton_improve_roots(prev_roots, f ,derivative_f, eps):
    # метод ньютона, улучщающий найденные корни
    roots = list()
    root_counts = [0] * len(prev_roots)

    for i, prev_root in enumerate(prev_roots):
        while True:
            root_counts[i] += 1

            x = prev_root - f(prev_root) / derivative_f(prev_root)

            if math.fabs(x - prev_root) < eps:
                break

            prev_root = x

        roots.append(x)

    return (roots, root_counts)

def newton(sections, f, derivative_f, eps):
    # метод ньютона, находящий корни из отрезков, сжатых бисекцией
    roots = list()
    root_counts = [0] * len(sections)

    for i, section in enumerate(sections):
        prev_x = section[1]

        while True:
            root_counts[i] += 1

            x = prev_x - f(prev_x) / derivative_f(prev_x)

            if math.fabs(x - prev_x) < eps:
                break

            prev_x = x

        roots.append(x)

    return (roots, root_counts)

if __name__ == "__main__":
    result_sections, section_counts = bisection(sections, f, BISECTION_SECTION_LENGTH)

    discrete_newton_roots, discrete_newton_root_counts = discrete_newton(result_sections, f, \
        DISCRETE_NEWTON_H, DISCRETE_NEWTON_EPS)

    newton_roots, newton_root_counts = newton(result_sections, f, derivative_f, NEWTON_EPS)
    newton_improve_roots, newton_improve_root_counts = newton_improve_roots(discrete_newton_roots, f, \
        derivative_f, NEWTON_EPS)

    with open("report.txt", "a") as report_file:
        report_file.write("\n")
        report_file.write("BISECTION (TASK 6):\n")
        report_file.write("Segments of separation: ")
        report_file.write(str(result_sections))
        report_file.write("\n")
        report_file.write("Number of steps for each segment: ")
        report_file.write(str(section_counts))
        report_file.write("\n")
        report_file.write("BISECTION_SECTION_LENGTH = " + str(BISECTION_SECTION_LENGTH))
        report_file.write("\n")
        report_file.write("\n")
        
        report_file.write("\n")
        report_file.write("DISCRETE NEWTON (TASK 7):\n")
        report_file.write("Roots: ")
        report_file.write(str(discrete_newton_roots))
        report_file.write("\n")
        report_file.write("Number of steps for each root: ")
        report_file.write(str(discrete_newton_root_counts))
        report_file.write("\n")
        report_file.write("EPS = " + str(DISCRETE_NEWTON_EPS) + ", H = " + str(DISCRETE_NEWTON_H))
        report_file.write("\n")
        report_file.write("\n")

        report_file.write("\n")
        report_file.write("IMPROVE ROOTS USING NEWTON (TASK 8):\n")
        report_file.write("Roots: ")
        report_file.write(str(newton_improve_roots))
        report_file.write("\n")
        report_file.write("Number of steps for each root: ")
        report_file.write(str(newton_improve_root_counts))
        report_file.write("\n")
        report_file.write("EPS = " + str(NEWTON_EPS))
        report_file.write("\n")
        report_file.write("\n")

        report_file.write("\n")
        report_file.write("NEWTON TO SEGMENTS OF SEPARATION (TASK 8):\n")
        report_file.write("Roots: ")
        report_file.write(str(newton_roots))
        report_file.write("\n")
        report_file.write("Number of steps for each root: ")
        report_file.write(str(newton_root_counts))
        report_file.write("\n")
        report_file.write("EPS = " + str(NEWTON_EPS))
        report_file.write("\n")

    draw_function(f, -26, 20, 0.01)


