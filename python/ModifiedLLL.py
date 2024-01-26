from decimal import *

def zero_judge(x):
    for i in x:
        if i != 0:
            return False
    return True

#ベクトルの和
def add(v1, v2):
    result = [0] * len(v1)
    for i in range(len(v1)):
        result[i] = Decimal(v1[i]) + Decimal(v2[i])
    return result

#ベクトルの差
def sub(v1, v2):
    result = [0] * len(v1) 
    for i in range(len(v1)):
        result[i] = Decimal(v1[i]) - Decimal(v2[i])
    return result

#ベクトルの積(定数倍)
def scalar(a, v1):
    result = [0] * len(v1)
    for i in range(len(v1)):
        result[i] = a * Decimal(v1[i]) 
    return result

#ベクトルの積(内積)
def inner_product(v1, v2):
    result = Decimal(0)
    for i in range(len(v1)):
        result += Decimal(v1[i]) * Decimal(v2[i])
    return result

#サイズ基底簡約
def sizeReduce(b,mu,i,j):
    if abs(mu[i][j]) > 0.5:
        q = round(mu[i][j])
        b[i] = sub(b[i], scalar(q,b[j]))
        for l in range(j+1):
            mu[i][l] -= q * Decimal(mu[j][l])
    return b,mu
    
def MLLL(b,delta):
    n = len(b[0]) 
    h = len(b) 
    b_gso = [([0] * n) for i in range(h)]
    Bi = [0] * h
    mu = [([0] * h) for i in range(h+1)]

    z = h - 1
    g = 0
    while g <= z:
        if zero_judge(b[g]):
            if g < z:
                b[g], b[z] = b[z], b[g]
            z -= 1

        b_gso[g] = b[g]
        mu[g][g] = 1
        for j in range(g):
            mu[g][j] = inner_product(b[g], b_gso[j]) / Bi[j]
            b_gso[g] = sub(b_gso[g], scalar(mu[g][j],b_gso[j]))
        Bi[g] = inner_product(b_gso[g], b_gso[g])

        if g == 0:
            g = 1 
        else:
            l = g
            k = g
            startagain = False
            while (k <= l) and (not startagain):
                b,mu = sizeReduce(b,mu,k,k-1)
                nu = mu[k][k-1]
                B = Bi[k] + nu**2 * Bi[k-1]

                if B >= Decimal(delta) * Bi[k-1]:
                    for j in reversed(range(k-1)):
                        b,mu = sizeReduce(b,mu,k,j)
                    k += 1
                else:
                    if zero_judge(b[k]):
                        if k < z:
                            b[k], b[z] = b[z], b[k]
                        z -= 1
                        g = k
                        startagain = True
                    else:
                        b[k-1], b[k] = b[k], b[k-1]
                        for j in range(k-1):
                            mu[k][j], mu[k-1][j] = mu[k-1][j], mu[k][j]
                        if B != 0:
                            if Bi[k] == 0:
                                Bi[k-1] = B
                                b_gso[k-1] = scalar(nu,b_gso[k-1])
                                mu[k][k-1] = 1 / nu
                                for i in range(k+1, l+1):
                                    mu[i][k-1] = mu[i][k-1] / nu
                            else:
                                t = Bi[k-1] / B
                                mu[k][k-1] = nu * t
                                w = b_gso[k-1]
                                b_gso[k-1] = add(b_gso[k], scalar(nu,w))
                                Bi[k-1] = B
                                if k <= l:
                                    b_gso[k] = sub(scalar((Bi[k]/B),w) , scalar(mu[k][k-1],b_gso[k]))
                                    Bi[k] = Bi[k] * t

                                for i in range(k+1, l+1):
                                    t = mu[i][k]
                                    mu[i][k] = mu[i][k-1] - nu * t
                                    mu[i][k-1] = t + mu[k][k-1] * mu[i][k]
                        else:
                            Bi[k], Bi[k-1] = Bi[k-1], Bi[k]
                            b_gso[k], b_gso[k-1] = b_gso[k-1], b_gso[k]
                            for i in range(k+1, l+1):
                                mu[i][k], mu[i][k-1] = mu[i][k-1], mu[i][k]
                        k = max(k-1,1)
            if not startagain:
                g += 1
    return [list(map(int, x)) for x in b]

b1 = [
      [388, 417, 417, -86],
      [-672, -73, -121, 944],
      [-689, 379, 724, 653],
      [-179, 96, -24, 978],
      [508, -705, 173, -343]]

b2 = [
      [-696, -186, 661, -727],
      [-760, -106, -775, 659],
      [552, 6, 9, 726],
      [-160, -439, -544, 365],
      [307, -526, 862, 396],
      [117, -94, 472, 138]]


result1 = MLLL(b1,delta=0.75)
result2 = MLLL(b2,delta=0.75)

for i in range(len(b1)):
    print(f"{result1[i]}")

print(f"\n")

for i in range(len(b2)):
    print(f"{result2[i]}")
