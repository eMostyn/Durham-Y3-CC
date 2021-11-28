import numpy as np
import math
import os
import random


def point_addition(P, Q, p, a):
    if P == [0, 0]:
        return Q
    if Q == [0, 0]:
        return P
    x1 = P[0]
    y1 = P[1]
    x2 = Q[0]
    y2 = Q[1]
    if x1 == x2 and y1 != y2:
        # p1 + -p1 == 0
        return [0, 0]
    if x1 == x2:
        # p1 + p1: use tangent line of p1 as (p1,p1) line
        l = (3 * x1 * x1 + a) * inv(2 * y1, p) % p
        #
    else:
        l = (y2 - y1) * inv(x2 - x1, p) % p
    x3 = ((l ** 2) - x1 - x2) % p
    y3 = (l * (x1 - x3) - y1) % p

    return [x3, y3]


def mult(p, n, P, a):
    r = (0, 0)
    if n == 2:
        return double(P, p, a)
    elif n > 2:
        multiplied = double(P, p, a)
        for i in range(n - 2):
            multiplied = point_addition(multiplied, P, p, a)
    return multiplied


def inv(n, q):
    for i in range(q):
        if (n * i) % q == 1:
            return i


def add(p, a, b):
    return (a + b) % p


def multiply(p, a, b):
    return (a * b) % p


def double(P, p, a):
    x1 = P[0]
    y1 = P[1]
    l = ((3 * (x1 ** 2) + a) * inv(2 * y1, p)) % p
    x3 = ((l) ** 2 - (2 * x1)) % p
    y3 = (l * (x1 - x3) - y1) % p
    return [x3, y3]


def generate_points(p, a, b):
    curve = [(p, p)]
    x_curve = []
    for x in range(p):
        x2 = (x ** 3 + a * x + b) % p
        x_curve.append((x, x2))
    y_curve = []
    for y in range(p):
        y2 = (y ** 2) % p
        y_curve.append((y, y2))
    y_curve.sort(key=lambda x: x[1])
    points = []
    for x in range(len(x_curve)):
        for y in range(0, len(y_curve)):
            if y_curve[y][1] > x_curve[x][1]:
                break
            if y_curve[y][1] == x_curve[x][1]:
                points.append((x_curve[x][0], y_curve[y][0]))
    return points


def H(x, points, partititon_1, partition_2):
    partition_1_num = points[partititon_1][0]
    partition_2_num = points[partition_2][0]
    if x < partition_1_num:
        return 0
    elif x < partition_2_num:
        return 1
    return 2


def F(X, x, As, Bs, P, p, a):
    temp = point_addition(X, (mult(p, As[x], P, a)), p, a)
    return point_addition(temp, mult(p, Bs[x], P, a), p, a)


def F2(X, x, P, p, a, c, d, n):
    if x == 0:
        return point_addition(X, P, p, a), add(n, c, 1), d
    elif x == 1:
        return mult(p, 2, X, a), multiply(n, 2, c), multiply(n, 2, d)
    else:
        return point_addition(X, Q, p, a), c, add(n, d, 1)


def pol_rho(p, a, b, P, n, Q):
    if os.path.exists("points.txt"):
        points = np.loadtxt("points.txt", dtype="int")
    else:
        points = generate_points(p, a, b)
        np.savetxt("points.txt", points, fmt="%d")
    num_points = len(points)
    partition_1 = num_points // 3
    partition_2 = (2 * num_points) // 3
    s_1 = points[:partition_1]
    s_2 = points[partition_1:partition_2]
    s_3 = points[partition_2:]
    # print(point_addition(mult(p, 753, P, a), mult(p, 3598, Q, a), p, a))

    As = [random.randint(0, n - 1) for x in range(3)]
    print(As)
    Bs = [random.randint(0, n - 1) for x in range(3)]
    print(Bs)

    # j = H(1452, points, partition_1, partition_2)
    # LHS = point_addition(mult(p, 3025, P, a), mult(p, 13260, Q, a), p)
    c = 5
    d = 9
    c_dash = 0
    d_dash = 0
    X = point_addition(mult(p, c, P, a), mult(p, d, Q, a), p, a)
    j = H(X[0], points, partition_1, partition_2)
    X_dash = F2(X, j, As, Bs, P, p, a)
    print(X, X_dash)

    while X != X_dash:
        j = H(X[0], points, partition_1, partition_2)
        X = F(X, j, As, Bs, P, p, a)
        c_new, d_new = As[j], Bs[j]
        c = add(p, c, c_new)
        d = add(p, d, d_new)
        j_x_dash = H(X_dash[0], points, partition_1, partition_2)
        c_dash_new, d_dash_new = As[j_x_dash], Bs[j_x_dash]
        X_dash = F(X_dash, j_x_dash, As, Bs, P, p, a)
        c_dash = add(p, c, c_dash_new)
        d_dash = add(p, d, d_dash_new)
        j_x_dash = H(X_dash[0], points, partition_1, partition_2)
        c_dash_new, d_dash_new = As[j_x_dash], Bs[j_x_dash]
        X_dash = F(X_dash, j_x_dash, As, Bs, P, p, a)
        c_dash = add(p, c, c_dash_new)
        d_dash = add(p, d, d_dash_new)
        print(X, X_dash)

        # C = 5 D=9
    return X, X_dash, c, d, c_dash, d_dash


def pol_rho2(p, a, b, P, n, Q):
    if os.path.exists("points.txt"):
        points = np.loadtxt("points.txt", dtype="int")
    else:
        points = generate_points(p, a, b)
        np.savetxt("points.txt", points, fmt="%d")
    num_points = len(points)
    partition_1 = num_points // 3
    partition_2 = (2 * num_points) // 3
    s_1 = points[:partition_1]
    s_2 = points[partition_1:partition_2]
    s_3 = points[partition_2:]
    # print(point_addition(mult(p, 753, P, a), mult(p, 3598, Q, a), p, a))

    As = [random.randint(0, n - 1) for x in range(3)]
    print(As)
    Bs = [random.randint(0, n - 1) for x in range(3)]
    print(Bs)

    # j = H(1452, points, partition_1, partition_2)
    # LHS = point_addition(mult(p, 3025, P, a), mult(p, 13260, Q, a), p)
    c = 22
    d = 10
    c_dash = 4
    d_dash = 8
    X = point_addition(mult(p, c, P, a), mult(p, d, Q, a), p, a)
    # X_dash = point_addition(mult(p, c_dash, P, a), mult(p, d_dash, Q, a), p, a)
    Xs = [[X, c, d]]
    print(X)
    found = False
    while not found:
        j = H(X[0], points, partition_1, partition_2)
        X, c, d = F2(X, j, P, p, a, c, d, n)
        unzipped = list(zip(*Xs))
        if X in unzipped[0]:
            index = unzipped[0].index(X)
            return Xs[index], X, c, d
        else:
            Xs.append([X, c, d])
        # print(X)
    # j_x_dash = H(X_dash[0], points, partition_1, partition_2)
    # X_dash, c_dash, d_dash = F2(X_dash, j_x_dash, P, p, a, c, d, n)
    # j_x_dash = H(X_dash[0], points, partition_1, partition_2)
    # X_dash, c_dash, d_dash = F2(X_dash, j_x_dash, P, p, a, c, d, n)


# return X, X_dash, c, d, c_dash, d_dash


def check_valid(P, Q, c, d, c_dash, d_dash, a, p):
    print("LHS: ", point_addition(mult(p, c, P, a), mult(p, d, Q, a), p, a))
    print("RHS: ", point_addition(mult(p, c_dash, P, a), mult(p, d_dash, Q, a), p, a))


text = np.genfromtxt("exampleInputRho.txt", dtype="str", comments="#", delimiter="\n")[
    1:
]
parameters = {}
for item in text:
    parameters[item[:1]] = item[4:]
p = int(parameters["p"])
a = int(parameters["a"])
b = int(parameters["b"])
P = parameters["P"].strip("()").split(", ")
P = [int(x) for x in P]
n = int(parameters["n"])
Q = parameters["Q"].strip("()").split(", ")
Q = [int(x) for x in Q]

print(pol_rho2(p, a, b, P, n, Q))
# check_valid(P, Q, 6150, 379, 4480, 3022, a, p)
# 1266 11846
# 5735 5226
