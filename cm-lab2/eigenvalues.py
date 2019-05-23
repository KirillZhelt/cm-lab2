from numpy import linalg as la

if __name__ == "__main__":
    with open("matrix.txt", "r") as matrix_file:
        matrix = [list(map(float, line.rstrip().split(" "))) for line in matrix_file.readlines()]

    with open("report.txt", "a") as report_file:
        report_file.write(str(la.eig(matrix)[0]))
        report_file.write("\n")
        report_file.write("\n")
