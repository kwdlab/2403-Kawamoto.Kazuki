import numpy as np

def sizeReduce(b,mu,i,j):
    if abs(mu[i,j]) > 0.5:
        q = round(mu[i,j])
        b[i] = b[i] - q * b[j]
        for l in range(j+1):
            mu[i,l] -= q * mu[j,l]
    return b,mu

def MLLL(b,delta):
    n = b.ncols()
    h = b.nrows()
    GSOVector = Matrix(QQ,h,n)
    mu = Matrix(QQ,h)
    Bi = vector(QQ, h)

    z = h-1
    g = 0
    while g <= z:
        if np.all(b[g] == 0):
            if g < z:
                v = b[g] 
                b[g] = b[z]
                b[z] = v
            z -= 1

        GSOVector[g] = b[g]
        mu[g,g] = 1
        for j in range(g):
            mu[g,j] = b[g].inner_product(GSOVector[j]) / GSOVector[j].norm()**2
            GSOVector[g] -= GSOVector[j] * mu[g,j]
        Bi[g] = GSOVector[g].norm() ** 2

        if g == 0:
            g = 1 
        else:
            l = g
            k = g
            startagain = False
            while k <= l and not startagain:
                b,mu = sizeReduce(b,mu,k,k-1)
                nu = mu[k,k-1]
                B = Bi[k] + nu * nu * Bi[k-1]

                if B >= delta * Bi[k-1]:
                    for j in reversed(range(k-1)):
                        b,mu = sizeReduce(b,mu,k,j)
                    k += 1
                else:
                    if np.all(b[k]==0):
                        if k < z:
                            v = b[k]
                            b[k] = b[z]
                            b[z] = v
                        z -= 1
                        g = k
                        startagain = True
                    else:
                        v = b[k]
                        b[k] = b[k-1]
                        b[k-1] = v
                        for j in range(k-1):
                           v =  mu[k,j]
                           mu[k,j] = mu[k-1,j] 
                           mu[k-1,j] = v

                        if B != 0:
                            if Bi[k] == 0:
                                Bi[k-1] = B
                                GSOVector[k-1] *= nu 
                                mu[k,k-1] = 1 / nu
                                for i in range(k+1,l+1):
                                    mu[i,k-1] = mu[i,k-1] / nu
                            else:
                                t = Bi[k-1] / B
                                mu[k,k-1] = nu * t
                                w = GSOVector[k-1]
                                GSOVector[k-1] = GSOVector[k] +  w * nu
                                Bi[k-1] = B
                                if k <= l:
                                    GSOVector[k] = (w * Bi[k]/B) -  GSOVector[k] * mu[k,k-1]
                                    Bi[k] = Bi[k] * t

                                for i in range(k+1,l+1):
                                    t = mu[i,k]
                                    mu[i,k] = mu[i,k-1] - nu * t
                                    mu[i,k-1] = t + mu[k,k-1] * mu[i,k]
                        else:
                            v = Bi[k] 
                            Bi[k] = Bi[k-1]
                            Bi[k-1] = v

                            v = GSOVector[k]  
                            GSOVector[k] = GSOVector[k-1]
                            GSOVector[k-1] = v

                            for i in range(k+1,l+1):
                                v = mu[i,k] 
                                mu[i,k] = mu[i,k-1]
                                mu[i,k-1] = v
                        k = max(k-1, 1)
            if not startagain:
                g += 1
    return b

b1 = Matrix([
    [388, 417, 417, -86],
    [-672, -73, -121, 944],
    [-689, 379, 724, 653],
    [-179, 96, -24, 978],
    [508, -705, 173, -343]
    ])

b2 = Matrix([
    [-696, -186, 661, -727],
    [-760, -106, -775, 659],
    [552, 6, 9, 726],
    [-160, -439, -544, 365],
    [307, -526, 862, 396],
    [117, -94, 472, 138]
])

result1 = MLLL(b1,delta=0.75)
result2 = MLLL(b2,delta=0.75)

print(result1)
print(f"\n")
print(result2)
