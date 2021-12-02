import numpy as np
import os
import random
import math
import time


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
        return [0, 0]
    if x1 == x2:
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
    else:
        return ValueError


def inv(n, p):
    # for i in range(q):
    #     if (n * i) % q == 1:
    #         return i

    return pow(n, -1, p)


def add(p, a, b):
    return (a + b) % p


def sub(p, a, b):
    return (b - a) % p


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


def H(x, order):
    partition_1_num = order // 3
    partition_2_num = (2 * order) // 3
    if x < partition_1_num:
        return 0
    elif x < partition_2_num:
        return 1
    return 2


def F2(X, x, P, Q, p, a, c, d, n):
    if x == 0:
        return point_addition(X, P, p, a), add(n, c, 1), d
    elif x == 1:
        return mult(p, 2, X, a), multiply(n, 2, c), multiply(n, 2, d)
    else:
        return point_addition(X, Q, p, a), c, add(n, d, 1)


def pol_rho2(p, a, b, P, n, Q):
    c = 2
    d = 87
    X = point_addition(mult(p, c, P, a), mult(p, d, Q, a), p, a)

    Xs = [[X, c, d]]
    found = False
    while not found:
        j = H(X[0], n)
        X, c, d = F2(X, j, P, Q, p, a, c, d, n)
        unzipped = list(zip(*Xs))
        if X in unzipped[0]:
            index = unzipped[0].index(X)
            found = True
            c, d, c_dash, d_dash = Xs[index][1], Xs[index][2], c, d
        else:
            Xs.append([X, c, d])
    print(c, d, c_dash, d_dash)
    return c, d, c_dash, d_dash


def get_l(c, d, c_dash, d_dash, n, p, P, Q):
    gcd = math.gcd(d_dash - d, n)
    if gcd == 1:
        l = sub(n, c, c_dash) * (inv(sub(n, d_dash, d), n)) % n
    elif gcd > 1:
        for i in range(1, gcd + 1):
            new_mod = int((n * i) / gcd)
            try:
                l = (
                    sub(new_mod, c, c_dash)
                    * (inv(sub(new_mod, d_dash, d), new_mod))
                    % new_mod
                )
            except:
                pass
    return l


def print_output(parameters, c, d, c_dash, d_dash):
    with open("OutputBasicRho.txt", "w") as outfile:
        outfile.write("Input:\n")
        for key in parameters.keys():
            line = key + " = " + parameters[key] + "\n"
            outfile.write(line)
        outfile.write("\nCollision:\n")
        outfile.write("c = " + str(c) + "\n")
        outfile.write("d = " + str(d) + "\n")
        outfile.write("c' = " + str(c_dash) + "\n")
        outfile.write("d' = " + str(d_dash) + "\n")


def check_valid(P, Q, c, d, c_dash, d_dash, a, p):
    print(P, Q)
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

start = time.time()
c, d, c_dash, d_dash = pol_rho2(p, a, b, P, n, Q)
print("Found collision: ", time.time() - start)

print(get_l(c, d, c_dash, d_dash, n, p, P, Q))
print("Found l: ", time.time() - start)
