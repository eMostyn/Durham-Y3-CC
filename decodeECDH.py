import numpy as np
import os
import random
import math
import time
from Crypto.Cipher import DES
import des


def point_addition(P, Q, p, a):
    if P ==  None:
        return Q
    if Q == None:
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


def mult2(p, n, P, a):
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

def bits(n):
    """
    Generates the binary digits of n, starting
    from the least significant bit.

    bits(151) -> 1, 1, 1, 0, 1, 0, 0, 1
    """
    while n:
        yield n & 1
        n >>= 1

def mult(p, n, P,a):
    """
    Returns the result of n * x, computed using
    the double and add algorithm.
    """
    result = None
    addend = P
    for bit in bits(n):
        if bit:
            result = point_addition(result,addend,p,a)

        addend = double(addend, p, a)

    return result

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
    c_dash = 0
    d_dash = 0
    X = point_addition(mult(p, c, P, a), mult(p, d, Q, a), p, a)
    j = H(X[0], n)
    X_dash, c_dash, d_dash = F2(X, j, P, Q, p, a, c, d, n)
    # X_dash = point_addition(mult(p, c_dash, P, a), mult(p, d_dash, Q, a), p, a)
    found = False
    while X != X_dash:
        j = H(X[0], n)
        X, c, d = F2(X, j, P, Q, p, a, c, d, n)
        j_dash = H(X_dash[0], n)
        X_dash, c_dash, d_dash = F2(X_dash, j_dash, P, Q, p, a, c_dash, d_dash, n)
        j_dash = H(X_dash[0], n)
        X_dash, c_dash, d_dash = F2(X_dash, j_dash, P, Q, p, a, c_dash, d_dash, n)
    print(c, d, c_dash, d_dash)
    return c, d, c_dash, d_dash


def get_l(c, d, c_dash, d_dash, n,Q,P):
    print(d_dash-d,n)
    gcd = math.gcd(d_dash - d, n)
    print("gcd:", gcd)
    if gcd == 1:
        l = sub(n, c, c_dash) * (inv(sub(n, d_dash, d), n)) % n
    elif gcd > 1:
        # new_mod = n//gcd
        # k = sub(new_mod, c, c_dash) * (inv(sub(new_mod, d_dash, d), new_mod)) % new_mod
        # print(k)
        # print(Q,mult(p,k,P,a))

        # for i in range(0, gcd):
        #     num = (n*i)//gcd
        #     result = mult(p,num,P,a)
        #     print(Q,result)
        #     if Q == result:
        #         print("works")
        #         return num

         for i in range(0, gcd + 1):
            new_mod = int((n * i) / gcd)
            try:
                l = (
                    sub(new_mod, c, c_dash)
                    * (inv(sub(new_mod, d_dash, d), new_mod))
                    % new_mod
                )
                print(l)
                print(Q,mult(p,13333,P,a))
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


def pad(text):
    n = len(text) % 8
    return text + (b'\x00' * n)


def bitstring_to_bytes(s):
    return int(s, 2).to_bytes((len(s) + 7) // 8, byteorder='big')

def decode():
    key = '6714934996831608'

    #key = '10111110110110011001010100011001011010000110101111000'
    cipher_text = "3da46f7b6fa82f53153908bdadcc742ac38e8691e5208aa4bf6be47240c71e75180b9d1030a00810"
    cipher_bytes = bytes.fromhex(cipher_text)
    #key = bytes(6714934996831608)
    #key = str.encode(key)
    #000101111101101100110010101000110010110100 0011010  1111000
    #Added parity bits
    #0001011011101100110011010101010000110010011010000011010011110001
    key = bitstring_to_bytes("0001011011101100110011010101010000110010011010000011010011110001")
    print(key)

    # key = bytes.fromhex(hex(6714934996831608))
    # print(key)
    des1 = DES.new(key, DES.MODE_ECB)
    decrypted = des1.decrypt(cipher_bytes)
    #decrypted[5:].decode("utf-8") 


    key0 = des.DesKey(key)  
    decrypted = (key0.decrypt(cipher_bytes, padding=True)) 
    print(decrypted.decode("utf-8"))

decode()


# p = 20376993552394903
# a = 10
# b = 1
# P = (1983, 6761152449250519)
# n = 1852453970120513
# QA = (18586784116581871, 12161036958498472)
# QB = (18432261261031243, 11140924411855488)


p =16001
a =1
b =5
P =(1300, 16000)
QA = (12025, 11184)
QB = (11537, 9748)
n = 16190
# c = 391784786894884
# d = 475799047835362
# c_dash = 142130452030523
# d_dash = 1795121125780147

# p = 16001
# a = 10
# b = 1
# P = (1654, 7208)
# n = 8026
# Q = (5000, 1283)
# start = time.time()
c, d, c_dash, d_dash = pol_rho2(p, a, b, P, n, QB)
# print("Found collision: ", time.time() - start)
print(c,d,c_dash,d_dash)
l = get_l(c, d, c_dash, d_dash, n,QB,P)
print("l",l)
# da=1682779984167835
#print(mult(p, l, QB, a))
# # print(mult(p,l,QB,a))
# db= 428971283427559
# # print(QB,mult(p,l,P,a))
# print(mult(p,db,QA,a))
# print(mult(p,db,mult(p, da, P, a),a))
# print("Found l: ", time.time() - start)


# 391784786894884 475799047835362 142130452030523 1795121125780147
# l = get_l(391784786894884,475799047835362,142130452030523,1795121125780147,1852453970120513,QB,P)
# print("L",l)
# # check_valid(
# #     P, QA, 391784786894884, 475799047835362, 142130452030523, 1795121125780147, a, p
# # )
# # check_valid(
# #     P, QB, 1374007329837300, 464020786357595, 1374236167093169, 1204491816720654, a, p
# # )

# # QA: l = 1682779984167835
# # 
# # 1374007329837300 464020786357595 1374236167093169 1204491816720654
# l = get_l(1374007329837300,464020786357595,1374236167093169,1204491816720654,1852453970120513,QA,P)
# print("L",l)
# QB: l = 428971283427559
# 


#Key = [6714934996831608, 12073846457401645]
#6714934996831608


# c, d, c_dash, d_dash = pol_rho2(719, 130, 565, [312, 90], 233, [475, 662])
# check_valid([312, 90], [475, 662], c, d, c_dash, d_dash, 130, 719)
# print(get_l(c, d, c_dash, d_dash, 233, 719, [312, 90], [475, 662]))

# p = 1009
# P = [909, 601]
# Q = [134, 52]
# n = 1007
# a = 250
# b = 844
# start = time.time()
# c, d, c_dash, d_dash = pol_rho2(p, a, b, P, n, Q)
# check_valid(P, Q, c, d, c_dash, d_dash, a, p)
# print("Found collision: ", time.time() - start)

# l = get_l(c, d, c_dash, d_dash, n)
# print(l)
# print("2aha3",Q,mult(p,l,P,a))
