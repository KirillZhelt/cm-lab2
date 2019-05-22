from numpy import linalg as la

if __name__ == "__main__":
    with open("matrix.txt", "r") as matrix_file:
        matrix = [list(map(float, line.rstrip().split(" "))) for line in matrix_file.readlines()]

    print(la.eig(matrix)[0])