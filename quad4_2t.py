from numpy import array, sqrt, zeros, ix_
from scipy.linalg import det, inv
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve

def quad4(xy, properties):
    E = properties["E"]
    ν = properties["nu"]
    bx = properties["bx"]
    by = properties["by"]
    t = properties["t"]

    Eσ = E / (1 - ν ** 2) * array(
        [
            [1, ν, 0],
            [ν, 1, 0],
            [0, 0, (1 - ν) / 2]
        ])

    x0 = xy[0, 0]
    x1 = xy[1, 0]
    x2 = xy[2, 0]
    x3 = xy[3, 0]

    y0 = xy[0, 1]
    y1 = xy[1, 1]
    y2 = xy[2, 1]
    y3 = xy[3, 1]

    ke = zeros((8, 8))
    fe = zeros((8, 1))

    # Primer punto de Gauss de la regla 2x2
    # xi = 1.0 / sqrt(3)
    # eta = -1.0 / sqrt(3)
    # wi = 1.0
    # wj = 1.0

    gauss_rule = [
        (-1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
        (1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
        (1.0 / sqrt(3), 1.0 / sqrt(3), 1.0, 1.0),
        (-1.0 / sqrt(3), 1.0 / sqrt(3), 1.0, 1.0),
    ]

    for xi, eta, wi, wj in gauss_rule:

        # print(f"xi = {xi} eta = {eta}")

        x = x0 * (1 - eta) * (1 - xi) / 4 + x1 * (1 - eta) * (xi + 1) / 4 + x2 * (eta + 1) * (xi + 1) / 4 + x3 * (
                    1 - xi) * (eta + 1) / 4
        y = y0 * (1 - eta) * (1 - xi) / 4 + y1 * (1 - eta) * (xi + 1) / 4 + y2 * (eta + 1) * (xi + 1) / 4 + y3 * (
                    1 - xi) * (eta + 1) / 4
        dx_dxi = -x0 * (1 - eta) / 4 + x1 * (1 - eta) / 4 + x2 * (eta + 1) / 4 - x3 * (eta + 1) / 4
        dx_deta = -x0 * (1 - xi) / 4 - x1 * (xi + 1) / 4 + x2 * (xi + 1) / 4 + x3 * (1 - xi) / 4
        dy_dxi = -y0 * (1 - eta) / 4 + y1 * (1 - eta) / 4 + y2 * (eta + 1) / 4 - y3 * (eta + 1) / 4
        dy_deta = -y0 * (1 - xi) / 4 - y1 * (xi + 1) / 4 + y2 * (xi + 1) / 4 + y3 * (1 - xi) / 4

        dN0_dxi = eta / 4. - 1 / 4.
        dN0_deta = xi / 4. - 1 / 4.
        dN1_dxi = 1 / 4. - eta / 4.
        dN1_deta = -xi / 4. - 1 / 4.
        dN2_dxi = eta / 4. + 1 / 4.
        dN2_deta = xi / 4. + 1 / 4.
        dN3_dxi = -eta / 4. - 1 / 4.
        dN3_deta = 1 / 4. - xi / 4.

        # print(f"x = {x} y = {y}")

        J = array([
            [dx_dxi, dx_deta],
            [dy_dxi, dy_deta]
        ]).T

        detJ = det(J)

        if detJ <= 0.:
            print(f"FATAL! detJ <= 0...")
            exit(-1)

        Jinv = inv(J)

        # print(f"J = {J}")
        # print(f"detJ = {detJ}")

        dN0_dxy = Jinv @ array([dN0_dxi, dN0_deta])
        dN1_dxy = Jinv @ array([dN1_dxi, dN1_deta])
        dN2_dxy = Jinv @ array([dN2_dxi, dN2_deta])
        dN3_dxy = Jinv @ array([dN3_dxi, dN3_deta])

        # ε = B ue
        B = lil_matrix(zeros((3, 8)))
        B[0, 0] = dN0_dxy[0]
        B[1, 1] = dN0_dxy[1]
        B[2, 0] = dN0_dxy[1]
        B[2, 1] = dN0_dxy[0]
        B[0, 2] = dN1_dxy[0]
        B[1, 3] = dN1_dxy[1]
        B[2, 2] = dN1_dxy[1]
        B[2, 3] = dN1_dxy[0]
        B[0, 4] = dN2_dxy[0]
        B[1, 5] = dN2_dxy[1]
        B[2, 4] = dN2_dxy[1]
        B[2, 5] = dN2_dxy[0]
        B[0, 6] = dN3_dxy[0]
        B[1, 7] = dN3_dxy[1]
        B[2, 6] = dN3_dxy[1]
        B[2, 7] = dN3_dxy[0]

        # print(f"B = {B}")

        ke += t * wi * wj * B.T @ Eσ @ B * detJ



    return ke, fe

def quad9(xy, properties):
    E = properties["E"]
    ν = properties["nu"]
    bx = properties["bx"]
    by = properties["by"]
    t = properties["t"]

    Ex = properties["Ex"]
    Ey = properties["Ey"]
    νxy = properties["nuxy"]
    νyx = properties["nuyx"]

    #Eσ = E / (1 - ν ** 2) * array(
    #    [
    #        [1, ν, 0],
    #        [ν, 1, 0],
    #        [0, 0, (1 - ν) / 2]
    #    ])

    Eσ = (1 / (1 - νxy * νyx)) * array(
        [
            [Ex, νxy * Ey, 0],
            [νxy * Ey, Ey, 0],
            [0, 0, Ex * Ey * (1 - νxy * νyx) / (Ex + Ey + 2 * Ey * νxy)]
        ])

    #N0 = (eta - 1) * (xi - 1) * xi * eta / 4
    #N1 = (eta + 1) * (xi - 1) * xi * eta / 4
    #N2 = (eta + 1) * (xi + 1) * xi * eta / 4
    #N3 = (eta - 1) * (xi + 1) * xi * eta / 4
    #N4 = (1 - eta**2) * (xi - 1) * xi /2
    #N5 = (1 - xi**2) * (eta + 1) * eta /2
    #N6 = (1 - eta**2) * (xi + 1) * xi /2
    #N7 = (1 - xi**2) * (eta - 1) * eta /2
    #N8 = (1- xi**2) * (1 - eta**2)

    x0 = xy[0, 0]
    x1 = xy[1, 0]
    x2 = xy[2, 0]
    x3 = xy[3, 0]
    x4 = xy[4, 0]
    x5 = xy[5, 0]
    x6 = xy[6, 0]
    x7 = xy[7, 0]
    x8 = xy[8, 0]

    y0 = xy[0, 1]
    y1 = xy[1, 1]
    y2 = xy[2, 1]
    y3 = xy[3, 1]
    y4 = xy[4, 1]
    y5 = xy[5, 1]
    y6 = xy[6, 1]
    y7 = xy[7, 1]
    y8 = xy[8, 1]

    ke = zeros((18, 18))
    fe = zeros((18, 1))

    # Primer punto de Gauss de la regla 2x2
    # xi = 1.0 / sqrt(3)
    # eta = -1.0 / sqrt(3)
    # wi = 1.0
    # wj = 1.0

    gauss_rule = [
        (-sqrt(3/5), -sqrt(3/5),        5/9,        5/9), #1
        ( sqrt(3/5), -sqrt(3/5),        5/9,        5/9), #3
        ( sqrt(3/5),  sqrt(3/5),        5/9,        5/9), #9
        (-sqrt(3/5),  sqrt(3/5),        5/9,        5/9), #7
        (       0.0, -sqrt(3/5), sqrt(40)/9, sqrt(40)/9), #2
        ( sqrt(3/5),        0.0, sqrt(40)/9, sqrt(40)/9), #6
        (       0.0,  sqrt(3/5), sqrt(40)/9, sqrt(40)/9), #8
        (-sqrt(3/5),        0.0, sqrt(40)/9, sqrt(40)/9), #4
        (       0.0,        0.0,        8/9,        8/9), #5
    ]

    for Ξ, ζ, wi, wj in gauss_rule:

        # print(f"xi = {xi} eta = {eta}")

        x = x0 * Ξ * ζ * (Ξ - 1) * (ζ - 1) / 4 + \
            x1 * Ξ * ζ * (Ξ + 1) * (ζ - 1) / 4 + \
            x2 * Ξ * ζ * (Ξ + 1) * (ζ + 1) / 4 + \
            x3 * Ξ * ζ * (Ξ - 1) * (ζ + 1) / 4 + \
            x4 * ζ * (1 - Ξ ** 2) * (ζ - 1) / 2 + \
            x5 * Ξ * (1 - ζ ** 2) * (Ξ + 1) / 2 + \
            x6 * ζ * (1 - Ξ ** 2) * (ζ + 1) / 2 + \
            x7 * Ξ * (1 - ζ ** 2) * (Ξ - 1) / 2 + \
            x8 * (1 - Ξ ** 2) * (1 - ζ ** 2)

        y = y0 * Ξ * ζ * (Ξ - 1) * (ζ - 1) / 4 + \
            y1 * Ξ * ζ * (Ξ + 1) * (ζ - 1) / 4 + \
            y2 * Ξ * ζ * (Ξ + 1) * (ζ + 1) / 4 + \
            y3 * Ξ * ζ * (Ξ - 1) * (ζ + 1) / 4 + \
            y4 * ζ * (1 - Ξ ** 2) * (ζ - 1) / 2 + \
            y5 * Ξ * (1 - ζ ** 2) * (Ξ + 1) / 2 + \
            y6 * ζ * (1 - Ξ ** 2) * (ζ + 1) / 2 + \
            y7 * Ξ * (1 - ζ ** 2) * (Ξ - 1) / 2 + \
            y8 * (1 - Ξ ** 2) * (1 - ζ ** 2)


        dN0_dΞ = Ξ * ζ * (ζ - 1) / 4 + ζ * (Ξ - 1) * (ζ - 1) / 4
        dN0_dζ = Ξ * ζ * (Ξ - 1) / 4 + Ξ * (Ξ - 1) * (ζ - 1) / 4
        dN1_dΞ = Ξ * ζ * (ζ - 1) / 4 + ζ * (Ξ + 1) * (ζ - 1) / 4
        dN1_dζ = Ξ * ζ * (Ξ + 1) / 4 + Ξ * (Ξ + 1) * (ζ - 1) / 4
        dN2_dΞ = Ξ * ζ * (ζ + 1) / 4 + ζ * (Ξ + 1) * (ζ + 1) / 4
        dN2_dζ = Ξ * ζ * (Ξ + 1) / 4 + Ξ * (Ξ + 1) * (ζ + 1) / 4
        dN3_dΞ = Ξ * ζ * (ζ + 1) / 4 + ζ * (Ξ - 1) * (ζ + 1) / 4
        dN3_dζ = Ξ * ζ * (Ξ - 1) / 4 + Ξ * (Ξ - 1) * (ζ + 1) / 4
        dN4_dΞ = -Ξ * ζ * (ζ - 1)
        dN4_dζ = ζ * (1 - Ξ ** 2) / 2 + (1 / 2 - Ξ ** 2 / 2) * (ζ - 1)
        dN5_dΞ = Ξ * (1 - ζ ** 2) / 2 + (1 - ζ ** 2) * (Ξ / 2 + 1 / 2)
        dN5_dζ = -Ξ * ζ * (Ξ + 1)
        dN6_dΞ = -Ξ * ζ * (ζ + 1)
        dN6_dζ = ζ * (1 - Ξ ** 2) / 2 + (1 - Ξ ** 2) * (ζ / 2 + 1 / 2)
        dN7_dΞ = Ξ * (1 - ζ ** 2) / 2 + (1 / 2 - ζ ** 2 / 2) * (Ξ - 1)
        dN7_dζ = -Ξ * ζ * (Ξ - 1)
        dN8_dΞ = -2 * Ξ * (1 - ζ ** 2)
        dN8_dζ = -2 * ζ * (1 - Ξ ** 2)

        dx_dΞ = x0 * Ξ * ζ * (ζ - 1) / 4 + \
                x0 * ζ * (Ξ - 1) * (ζ - 1) / 4 + \
                x1 * Ξ * ζ * (ζ - 1) / 4 + \
                x1 * ζ * (Ξ + 1) * (ζ - 1) / 4 + \
                x2 * Ξ * ζ * (ζ + 1) / 4 + \
                x2 * ζ * (Ξ + 1) * (ζ + 1) / 4 + \
                x3 * Ξ * ζ * (ζ + 1) / 4 + \
                x3 * ζ * (Ξ - 1) * (ζ + 1) / 4 - \
                x4 * Ξ * ζ * (ζ - 1) + \
                x5 * Ξ * (1 - ζ ** 2) / 2 + \
                x5 * (1 - ζ ** 2) * (Ξ + 1) / 2 - \
                x6 * Ξ * ζ * (ζ + 1) + \
                x7 * Ξ * (1 - ζ ** 2) / 2 + \
                x7 * (1 - ζ ** 2) * (Ξ - 1) / 2 - \
                2 * x8 * Ξ * (1 - ζ ** 2)
        dx_dζ = x0 * Ξ * ζ * (Ξ - 1) / 4 + \
                x0 * Ξ * (Ξ - 1) * (ζ - 1) / 4 + \
                x1 * Ξ * ζ * (Ξ + 1) / 4 + \
                x1 * Ξ * (Ξ + 1) * (ζ - 1) / 4 + \
                x2 * Ξ * ζ * (Ξ + 1) / 4 + \
                x2 * Ξ * (Ξ + 1) * (ζ + 1) / 4 + \
                x3 * Ξ * ζ * (Ξ - 1) / 4 + \
                x3 * Ξ * (Ξ - 1) * (ζ + 1) / 4 + \
                x4 * ζ * (1 - Ξ ** 2) / 2 + \
                x4 * (1 - Ξ ** 2) * (ζ - 1) / 2 - \
                x5 * Ξ * ζ * (Ξ + 1) + \
                x6 * ζ * (1 - Ξ ** 2) / 2 + \
                x6 * (1 - Ξ ** 2) * (ζ + 1) / 2 - \
                x7 * Ξ * ζ * (Ξ - 1) - \
                2 * x8 * ζ * (1 - Ξ ** 2)
        dy_dΞ = y0 * Ξ * ζ * (ζ - 1) / 4 + \
                y0 * ζ * (Ξ - 1) * (ζ - 1) / 4 + \
                y1 * Ξ * ζ * (ζ - 1) / 4 + \
                y1 * ζ * (Ξ + 1) * (ζ - 1) / 4 + \
                y2 * Ξ * ζ * (ζ + 1) / 4 + \
                y2 * ζ * (Ξ + 1) * (ζ + 1) / 4 + \
                y3 * Ξ * ζ * (ζ + 1) / 4 + \
                y3 * ζ * (Ξ - 1) * (ζ + 1) / 4 - \
                y4 * Ξ * ζ * (ζ - 1) + \
                y5 * Ξ * (1 - ζ ** 2) / 2 + \
                y5 * (1 - ζ ** 2) * (Ξ + 1) / 2 - \
                y6 * Ξ * ζ * (ζ + 1) + \
                y7 * Ξ * (1 - ζ ** 2) / 2 + \
                y7 * (1 - ζ ** 2) * (Ξ - 1) / 2 - \
                2 * y8 * Ξ * (1 - ζ ** 2)
        dy_dζ = y0 * Ξ * ζ * (Ξ - 1) / 4 + \
                y0 * Ξ * (Ξ - 1) * (ζ - 1) / 4 + \
                y1 * Ξ * ζ * (Ξ + 1) / 4 + \
                y1 * Ξ * (Ξ + 1) * (ζ - 1) / 4 + \
                y2 * Ξ * ζ * (Ξ + 1) / 4 + \
                y2 * Ξ * (Ξ + 1) * (ζ + 1) / 4 + \
                y3 * Ξ * ζ * (Ξ - 1) / 4 + \
                y3 * Ξ * (Ξ - 1) * (ζ + 1) / 4 + \
                y4 * ζ * (1 - Ξ ** 2) / 2 + \
                y4 * (1 - Ξ ** 2) * (ζ - 1) / 2 - \
                y5 * Ξ * ζ * (Ξ + 1) + \
                y6 * ζ * (1 - Ξ ** 2) / 2 + \
                y6 * (1 - Ξ ** 2) * (ζ + 1) / 2 - \
                y7 * Ξ * ζ * (Ξ - 1) - \
                2 * y8 * ζ * (1 - Ξ ** 2)

        #dN0_dxi = eta * xi * (eta - 1) / 4 + eta * (eta - 1) * (xi - 1) / 4
        #dN0_deta = eta * xi * (xi - 1) / 4 + xi * (eta - 1) * (xi - 1) / 4
        #dN1_dxi = eta * xi * (eta + 1) / 4 + eta * (eta + 1) * (xi - 1) / 4
        #dN1_deta = eta * xi * (xi - 1) / 4 + xi * (eta + 1) * (xi - 1) / 4
        #dN2_dxi = eta * xi * (eta + 1) / 4 + eta * (eta + 1) * (xi + 1) / 4
        #dN2_deta = eta * xi * (xi + 1) / 4 + xi * (eta + 1) * (xi + 1) / 4
        #dN3_dxi = eta * xi * (eta - 1) / 4 + eta * (eta - 1) * (xi + 1) / 4
        #dN3_deta = eta * xi * (xi + 1) / 4 + xi * (eta - 1) * (xi + 1) / 4
        #dN4_dxi = xi * (1 - eta ** 2) / 2 + (1 / 2 - eta ** 2 / 2) * (xi - 1)
        #dN4_deta = -eta * xi * (xi - 1)
        #dN5_dxi = -eta * xi * (eta + 1)
        #dN5_deta = eta * (1 - xi ** 2) / 2 + (1 - xi ** 2) * (eta / 2 + 1 / 2)
        #dN6_dxi = xi * (1 - eta ** 2) / 2 + (1 - eta ** 2) * (xi / 2 + 1 / 2)
        #dN6_deta = -eta * xi * (xi + 1)
        #dN7_dxi = -eta * xi * (eta - 1)
        #dN7_deta = eta * (1 - xi ** 2) / 2 + (1 / 2 - xi ** 2 / 2) * (eta - 1)
        #dN8_dxi = -2 * xi * (1 - eta ** 2)
        #dN8_deta = -2 * eta * (1 - xi ** 2)

        # print(f"x = {x} y = {y}")

        J = array([
            [dx_dΞ, dx_dζ],
            [dy_dΞ, dy_dζ]
        ]).T

        detJ = det(J)
        #xy = xy.todense()
        #print(f"xy = {xy}")
        #print(f"J = {J}")
        #print(f"detJ = {detJ}")

        if detJ <= 0.:
            print(f"FATAL! detJ <= 0...")
            exit(-1)

        Jinv = inv(J)

        dN0_dxy = Jinv @ array([dN0_dΞ, dN0_dζ])
        dN1_dxy = Jinv @ array([dN1_dΞ, dN1_dζ])
        dN2_dxy = Jinv @ array([dN2_dΞ, dN2_dζ])
        dN3_dxy = Jinv @ array([dN3_dΞ, dN3_dζ])
        dN4_dxy = Jinv @ array([dN4_dΞ, dN4_dζ])
        dN5_dxy = Jinv @ array([dN5_dΞ, dN5_dζ])
        dN6_dxy = Jinv @ array([dN6_dΞ, dN6_dζ])
        dN7_dxy = Jinv @ array([dN7_dΞ, dN7_dζ])
        dN8_dxy = Jinv @ array([dN8_dΞ, dN8_dζ])

        # ε = B ue
        B = lil_matrix(zeros((3, 18)))
        B[0, 0] = dN0_dxy[0]
        B[1, 1] = dN0_dxy[1]
        B[2, 0] = dN0_dxy[1]
        B[2, 1] = dN0_dxy[0]
        B[0, 2] = dN1_dxy[0]
        B[1, 3] = dN1_dxy[1]
        B[2, 2] = dN1_dxy[1]
        B[2, 3] = dN1_dxy[0]
        B[0, 4] = dN2_dxy[0]
        B[1, 5] = dN2_dxy[1]
        B[2, 4] = dN2_dxy[1]
        B[2, 5] = dN2_dxy[0]
        B[0, 6] = dN3_dxy[0]
        B[1, 7] = dN3_dxy[1]
        B[2, 6] = dN3_dxy[1]
        B[2, 7] = dN3_dxy[0]
        B[0, 8] = dN4_dxy[0]
        B[1, 9] = dN4_dxy[1]
        B[2, 8] = dN4_dxy[1]
        B[2, 9] = dN4_dxy[0]
        B[0, 10] = dN5_dxy[0]
        B[1, 11] = dN5_dxy[1]
        B[2, 10] = dN5_dxy[1]
        B[2, 11] = dN5_dxy[0]
        B[0, 12] = dN6_dxy[0]
        B[1, 13] = dN6_dxy[1]
        B[2, 12] = dN6_dxy[1]
        B[2, 13] = dN6_dxy[0]
        B[0, 14] = dN7_dxy[0]
        B[1, 15] = dN7_dxy[1]
        B[2, 14] = dN7_dxy[1]
        B[2, 15] = dN7_dxy[0]
        B[0, 16] = dN8_dxy[0]
        B[1, 17] = dN8_dxy[1]
        B[2, 16] = dN8_dxy[1]
        B[2, 17] = dN8_dxy[0]

        # print(f"B = {B}")

        ke += t * wi * wj * B.T @ Eσ @ B * detJ



    return ke, fe


def quad4_post(xy, u_e, properties):
    E = properties["E"]
    ν = properties["nu"]
    bx = properties["bx"]
    by = properties["by"]
    t = properties["t"]

    # Podemos pasarle otros valores de xi y eta
    if "xi" in properties:
        xi = properties["xi"]
    else:
        xi = 0.0

    if "eta" in properties:
        eta = properties["eta"]
    else:
        eta = 0.0

    Eσ = E / (1 - ν ** 2) * array(
        [
            [1, ν, 0],
            [ν, 1, 0],
            [0, 0, (1 - ν) / 2]
        ])


    x0 = xy[0, 0]
    x1 = xy[1, 0]
    x2 = xy[2, 0]
    x3 = xy[3, 0]

    y0 = xy[0, 1]
    y1 = xy[1, 1]
    y2 = xy[2, 1]
    y3 = xy[3, 1]

    x = x0 * (1 - eta) * (1 - xi) / 4 + x1 * (1 - eta) * (xi + 1) / 4 + x2 * (eta + 1) * (xi + 1) / 4 + x3 * (
                1 - xi) * (eta + 1) / 4
    y = y0 * (1 - eta) * (1 - xi) / 4 + y1 * (1 - eta) * (xi + 1) / 4 + y2 * (eta + 1) * (xi + 1) / 4 + y3 * (
                1 - xi) * (eta + 1) / 4
    dx_dxi = -x0 * (1 - eta) / 4 + x1 * (1 - eta) / 4 + x2 * (eta + 1) / 4 - x3 * (eta + 1) / 4
    dx_deta = -x0 * (1 - xi) / 4 - x1 * (xi + 1) / 4 + x2 * (xi + 1) / 4 + x3 * (1 - xi) / 4
    dy_dxi = -y0 * (1 - eta) / 4 + y1 * (1 - eta) / 4 + y2 * (eta + 1) / 4 - y3 * (eta + 1) / 4
    dy_deta = -y0 * (1 - xi) / 4 - y1 * (xi + 1) / 4 + y2 * (xi + 1) / 4 + y3 * (1 - xi) / 4

    dN0_dxi = eta / 4. - 1 / 4.
    dN0_deta = xi / 4. - 1 / 4.
    dN1_dxi = 1 / 4. - eta / 4.
    dN1_deta = -xi / 4. - 1 / 4.
    dN2_dxi = eta / 4. + 1 / 4.
    dN2_deta = xi / 4. + 1 / 4.
    dN3_dxi = -eta / 4. - 1 / 4.
    dN3_deta = 1 / 4. - xi / 4.

    # print(f"x = {x} y = {y}")

    J = array([
        [dx_dxi, dx_deta],
        [dy_dxi, dy_deta]
    ]).T

    detJ = det(J)

    if detJ <= 0.:
        print(f"FATAL! detJ <= 0...")
        exit(-1)

    Jinv = inv(J)

    # print(f"J = {J}")
    # print(f"detJ = {detJ}")

    dN0_dxy = Jinv @ array([dN0_dxi, dN0_deta])
    dN1_dxy = Jinv @ array([dN1_dxi, dN1_deta])
    dN2_dxy = Jinv @ array([dN2_dxi, dN2_deta])
    dN3_dxy = Jinv @ array([dN3_dxi, dN3_deta])

    # ε = B ue
    B = zeros((3, 8))
    B[0, 0] = dN0_dxy[0]
    B[1, 1] = dN0_dxy[1]
    B[2, 0] = dN0_dxy[1]
    B[2, 1] = dN0_dxy[0]
    B[0, 2] = dN1_dxy[0]
    B[1, 3] = dN1_dxy[1]
    B[2, 2] = dN1_dxy[1]
    B[2, 3] = dN1_dxy[0]
    B[0, 4] = dN2_dxy[0]
    B[1, 5] = dN2_dxy[1]
    B[2, 4] = dN2_dxy[1]
    B[2, 5] = dN2_dxy[0]
    B[0, 6] = dN3_dxy[0]
    B[1, 7] = dN3_dxy[1]
    B[2, 6] = dN3_dxy[1]
    B[2, 7] = dN3_dxy[0]

    u_e = lil_matrix(u_e)
    B = lil_matrix(B)

    #print(B)
    #print(u_e)

    B = B.todense()
    u_e = u_e.toarray()

    #print(B)
    #print(u_e)

    ε = B @ u_e.T
    σ = Eσ @ ε

    return ε, σ

def quad4n_post(xy, u_e, properties):
    E = properties["E"]
    ν = properties["nu"]
    bx = properties["bx"]
    by = properties["by"]
    t = properties["t"]

    # Podemos pasarle otros valores de xi y eta
    if "xi" in properties:
        xi = properties["xi"]
    else:
        xi = 0.0

    if "eta" in properties:
        eta = properties["eta"]
    else:
        eta = 0.0

    Eσ = E / (1 - ν ** 2) * array(
        [
            [1, ν, 0],
            [ν, 1, 0],
            [0, 0, (1 - ν) / 2]
        ])


    x0 = xy[0, 0]
    x1 = xy[1, 0]
    x2 = xy[2, 0]
    x3 = xy[3, 0]

    y0 = xy[0, 1]
    y1 = xy[1, 1]
    y2 = xy[2, 1]
    y3 = xy[3, 1]

    x = x0 * (1 - eta) * (1 - xi) / 4 + x1 * (1 - eta) * (xi + 1) / 4 + x2 * (eta + 1) * (xi + 1) / 4 + x3 * (
                1 - xi) * (eta + 1) / 4
    y = y0 * (1 - eta) * (1 - xi) / 4 + y1 * (1 - eta) * (xi + 1) / 4 + y2 * (eta + 1) * (xi + 1) / 4 + y3 * (
                1 - xi) * (eta + 1) / 4
    dx_dxi = -x0 * (1 - eta) / 4 + x1 * (1 - eta) / 4 + x2 * (eta + 1) / 4 - x3 * (eta + 1) / 4
    dx_deta = -x0 * (1 - xi) / 4 - x1 * (xi + 1) / 4 + x2 * (xi + 1) / 4 + x3 * (1 - xi) / 4
    dy_dxi = -y0 * (1 - eta) / 4 + y1 * (1 - eta) / 4 + y2 * (eta + 1) / 4 - y3 * (eta + 1) / 4
    dy_deta = -y0 * (1 - xi) / 4 - y1 * (xi + 1) / 4 + y2 * (xi + 1) / 4 + y3 * (1 - xi) / 4

    dN0_dxi = eta / 4. - 1 / 4.
    dN0_deta = xi / 4. - 1 / 4.
    dN1_dxi = 1 / 4. - eta / 4.
    dN1_deta = -xi / 4. - 1 / 4.
    dN2_dxi = eta / 4. + 1 / 4.
    dN2_deta = xi / 4. + 1 / 4.
    dN3_dxi = -eta / 4. - 1 / 4.
    dN3_deta = 1 / 4. - xi / 4.

    # print(f"x = {x} y = {y}")

    J = array([
        [dx_dxi, dx_deta],
        [dy_dxi, dy_deta]
    ]).T

    detJ = det(J)

    if detJ <= 0.:
        print(f"FATAL! detJ <= 0...")
        exit(-1)

    Jinv = inv(J)

    # print(f"J = {J}")
    # print(f"detJ = {detJ}")

    dN0_dxy = Jinv @ array([dN0_dxi, dN0_deta])
    dN1_dxy = Jinv @ array([dN1_dxi, dN1_deta])
    dN2_dxy = Jinv @ array([dN2_dxi, dN2_deta])
    dN3_dxy = Jinv @ array([dN3_dxi, dN3_deta])

    # ε = B ue
    B = zeros((3, 8))
    B[0, 0] = dN0_dxy[0]
    B[1, 1] = dN0_dxy[1]
    B[2, 0] = dN0_dxy[1]
    B[2, 1] = dN0_dxy[0]
    B[0, 2] = dN1_dxy[0]
    B[1, 3] = dN1_dxy[1]
    B[2, 2] = dN1_dxy[1]
    B[2, 3] = dN1_dxy[0]
    B[0, 4] = dN2_dxy[0]
    B[1, 5] = dN2_dxy[1]
    B[2, 4] = dN2_dxy[1]
    B[2, 5] = dN2_dxy[0]
    B[0, 6] = dN3_dxy[0]
    B[1, 7] = dN3_dxy[1]
    B[2, 6] = dN3_dxy[1]
    B[2, 7] = dN3_dxy[0]

    u_e = lil_matrix(u_e)
    B = lil_matrix(B)

    #print(B)
    #print(u_e)

    B = B.todense()
    u_e = u_e.toarray()

    #print(B)
    #print(u_e)

    εn = B @ u_e.T
    σn = Eσ @ ε

    return εn, σn


def quad9n_post(xy, u_e, properties):
    E = properties["E"]
    ν = properties["nu"]
    bx = properties["bx"]
    by = properties["by"]
    t = properties["t"]

    Ex = properties["Ex"]
    Ey = properties["Ey"]
    νxy = properties["nuxy"]
    νyx = properties["nuyx"]

    # Podemos pasarle otros valores de xi y eta
    if "xi" in properties:
        Ξ = properties["xi"]
    else:
        Ξ = 0.0

    if "eta" in properties:
        ζ = properties["eta"]
    else:
        ζ = 0.0

    #Eσ = E / (1 - ν ** 2) * array(
    #    [
    #        [1, ν, 0],
    #        [ν, 1, 0],
    #        [0, 0, (1 - ν) / 2]
    #    ])
#############
    Eσ = (1/(1 - νxy*νyx))*array(
        [
            [Ex, νxy*Ey, 0],
            [νxy*Ey, Ey, 0],
            [0, 0, Ex*Ey*(1 - νxy*νyx) / (Ex + Ey + 2*Ey*νxy)]
        ])
###############
    x0 = xy[0, 0]
    x1 = xy[1, 0]
    x2 = xy[2, 0]
    x3 = xy[3, 0]
    x4 = xy[4, 0]
    x5 = xy[5, 0]
    x6 = xy[6, 0]
    x7 = xy[7, 0]
    x8 = xy[8, 0]

    y0 = xy[0, 1]
    y1 = xy[1, 1]
    y2 = xy[2, 1]
    y3 = xy[3, 1]
    y4 = xy[4, 1]
    y5 = xy[5, 1]
    y6 = xy[6, 1]
    y7 = xy[7, 1]
    y8 = xy[8, 1]

    x = x0 * Ξ * ζ * (Ξ - 1) * (ζ - 1) / 4 + \
        x1 * Ξ * ζ * (Ξ + 1) * (ζ - 1) / 4 + \
        x2 * Ξ * ζ * (Ξ + 1) * (ζ + 1) / 4 + \
        x3 * Ξ * ζ * (Ξ - 1) * (ζ + 1) / 4 + \
        x4 * ζ * (1 - Ξ ** 2) * (ζ - 1) / 2 + \
        x5 * Ξ * (1 - ζ ** 2) * (Ξ + 1) / 2 + \
        x6 * ζ * (1 - Ξ ** 2) * (ζ + 1) / 2 + \
        x7 * Ξ * (1 - ζ ** 2) * (Ξ - 1) / 2 + \
        x8 * (1 - Ξ ** 2) * (1 - ζ ** 2)

    y = y0 * Ξ * ζ * (Ξ - 1) * (ζ - 1) / 4 + \
        y1 * Ξ * ζ * (Ξ + 1) * (ζ - 1) / 4 + \
        y2 * Ξ * ζ * (Ξ + 1) * (ζ + 1) / 4 + \
        y3 * Ξ * ζ * (Ξ - 1) * (ζ + 1) / 4 + \
        y4 * ζ * (1 - Ξ ** 2) * (ζ - 1) / 2 + \
        y5 * Ξ * (1 - ζ ** 2) * (Ξ + 1) / 2 + \
        y6 * ζ * (1 - Ξ ** 2) * (ζ + 1) / 2 + \
        y7 * Ξ * (1 - ζ ** 2) * (Ξ - 1) / 2 + \
        y8 * (1 - Ξ ** 2) * (1 - ζ ** 2)

    dx_dΞ = x0 * Ξ * ζ * (ζ - 1) / 4 + \
            x0 * ζ * (Ξ - 1) * (ζ - 1) / 4 + \
            x1 * Ξ * ζ * (ζ - 1) / 4 + \
            x1 * ζ * (Ξ + 1) * (ζ - 1) / 4 + \
            x2 * Ξ * ζ * (ζ + 1) / 4 + \
            x2 * ζ * (Ξ + 1) * (ζ + 1) / 4 + \
            x3 * Ξ * ζ * (ζ + 1) / 4 + \
            x3 * ζ * (Ξ - 1) * (ζ + 1) / 4 - \
            x4 * Ξ * ζ * (ζ - 1) + \
            x5 * Ξ * (1 - ζ ** 2) / 2 + \
            x5 * (1 - ζ ** 2) * (Ξ + 1) / 2 - \
            x6 * Ξ * ζ * (ζ + 1) + \
            x7 * Ξ * (1 - ζ ** 2) / 2 + \
            x7 * (1 - ζ ** 2) * (Ξ - 1) / 2 - \
            2 * x8 * Ξ * (1 - ζ ** 2)
    dx_dζ = x0 * Ξ * ζ * (Ξ - 1) / 4 + \
            x0 * Ξ * (Ξ - 1) * (ζ - 1) / 4 + \
            x1 * Ξ * ζ * (Ξ + 1) / 4 + \
            x1 * Ξ * (Ξ + 1) * (ζ - 1) / 4 + \
            x2 * Ξ * ζ * (Ξ + 1) / 4 + \
            x2 * Ξ * (Ξ + 1) * (ζ + 1) / 4 + \
            x3 * Ξ * ζ * (Ξ - 1) / 4 + \
            x3 * Ξ * (Ξ - 1) * (ζ + 1) / 4 + \
            x4 * ζ * (1 - Ξ ** 2) / 2 + \
            x4 * (1 - Ξ ** 2) * ( ζ - 1) / 2 - \
            x5 * Ξ * ζ * (Ξ + 1) + \
            x6 * ζ * (1 - Ξ ** 2) / 2 + \
            x6 * (1 - Ξ ** 2) * ( ζ + 1) / 2 - \
            x7 * Ξ * ζ * (Ξ - 1) - \
            2 * x8 * ζ * (1 - Ξ ** 2)
    dy_dΞ = y0 * Ξ * ζ * (ζ - 1) / 4 + \
            y0 * ζ * (Ξ - 1) * (ζ - 1) / 4 + \
            y1 * Ξ * ζ * (ζ - 1) / 4 + \
            y1 * ζ * (Ξ + 1) * ( ζ - 1) / 4 + \
            y2 * Ξ * ζ * (ζ + 1) / 4 + \
            y2 * ζ * (Ξ + 1) * (ζ + 1) / 4 + \
            y3 * Ξ * ζ * (ζ + 1) / 4 + \
            y3 * ζ * (Ξ - 1) * (ζ + 1) / 4 - \
            y4 * Ξ * ζ * (ζ - 1) + \
            y5 * Ξ * (1 - ζ ** 2) / 2 + \
            y5 * (1 - ζ ** 2) * (Ξ + 1) / 2 - \
            y6 * Ξ * ζ * (ζ + 1) + \
            y7 * Ξ * (1 - ζ ** 2) / 2 + \
            y7 * (1 - ζ ** 2) * (Ξ - 1) / 2 - \
            2 * y8 * Ξ * (1 - ζ ** 2)
    dy_dζ = y0 * Ξ * ζ * (Ξ - 1) / 4 + \
            y0 * Ξ * (Ξ - 1) * (ζ - 1) / 4 + \
            y1 * Ξ * ζ * (Ξ + 1) / 4 + \
            y1 * Ξ * (Ξ + 1) * (ζ - 1) / 4 + \
            y2 * Ξ * ζ * (Ξ + 1) / 4 + \
            y2 * Ξ * (Ξ + 1) * (ζ + 1) / 4 + \
            y3 * Ξ * ζ * (Ξ - 1) / 4 + \
            y3 * Ξ * (Ξ - 1) * (ζ + 1) / 4 + \
            y4 * ζ * (1 - Ξ ** 2) / 2 + \
            y4 * (1 - Ξ ** 2) * (ζ - 1) / 2 - \
            y5 * Ξ * ζ * (Ξ + 1) + \
            y6 * ζ * (1 - Ξ ** 2) / 2 + \
            y6 * (1 - Ξ ** 2) * (ζ + 1) / 2 - \
            y7 * Ξ * ζ * (Ξ - 1) - \
            2 * y8 * ζ * (1 - Ξ ** 2)

    dN0_dΞ = Ξ * ζ * (ζ - 1) / 4 + ζ * (Ξ - 1) * (ζ - 1) / 4
    dN0_dζ = Ξ * ζ * (Ξ - 1) / 4 + Ξ * (Ξ - 1) * (ζ - 1) / 4
    dN1_dΞ = Ξ * ζ * (ζ - 1) / 4 + ζ * (Ξ + 1) * (ζ - 1) / 4
    dN1_dζ = Ξ * ζ * (Ξ + 1) / 4 + Ξ * (Ξ + 1) * (ζ - 1) / 4
    dN2_dΞ = Ξ * ζ * (ζ + 1) / 4 + ζ * (Ξ + 1) * (ζ + 1) / 4
    dN2_dζ = Ξ * ζ * (Ξ + 1) / 4 + Ξ * (Ξ + 1) * (ζ + 1) / 4
    dN3_dΞ = Ξ * ζ * (ζ + 1) / 4 + ζ * (Ξ - 1) * (ζ + 1) / 4
    dN3_dζ = Ξ * ζ * (Ξ - 1) / 4 + Ξ * (Ξ - 1) * (ζ + 1) / 4
    dN4_dΞ = -Ξ * ζ * (ζ - 1)
    dN4_dζ = ζ * (1 - Ξ ** 2) / 2 + (1 / 2 - Ξ ** 2 / 2) * (ζ - 1)
    dN5_dΞ = Ξ * (1 - ζ ** 2) / 2 + (1 - ζ ** 2) * (Ξ / 2 + 1 / 2)
    dN5_dζ = -Ξ * ζ * (Ξ + 1)
    dN6_dΞ = -Ξ * ζ * (ζ + 1)
    dN6_dζ = ζ * (1 - Ξ ** 2) / 2 + (1 - Ξ ** 2) * (ζ / 2 + 1 / 2)
    dN7_dΞ = Ξ * (1 - ζ ** 2) / 2 + (1 / 2 - ζ ** 2 / 2) * (Ξ - 1)
    dN7_dζ = -Ξ * ζ * (Ξ - 1)
    dN8_dΞ = -2 * Ξ * (1 - ζ ** 2)
    dN8_dζ = -2 * ζ * (1 - Ξ ** 2)

    # print(f"x = {x} y = {y}")

    J = array([
        [dx_dΞ, dx_dζ],
        [dy_dΞ, dy_dζ]
    ]).T

    detJ = det(J)

    if detJ <= 0.:
        print(f"FATAL! detJ <= 0...")
        exit(-1)

    Jinv = inv(J)

    # print(f"J = {J}")
    # print(f"detJ = {detJ}")

    dN0_dxy = Jinv @ array([dN0_dΞ, dN0_dζ])
    dN1_dxy = Jinv @ array([dN1_dΞ, dN1_dζ])
    dN2_dxy = Jinv @ array([dN2_dΞ, dN2_dζ])
    dN3_dxy = Jinv @ array([dN3_dΞ, dN3_dζ])
    dN4_dxy = Jinv @ array([dN4_dΞ, dN4_dζ])
    dN5_dxy = Jinv @ array([dN5_dΞ, dN5_dζ])
    dN6_dxy = Jinv @ array([dN6_dΞ, dN6_dζ])
    dN7_dxy = Jinv @ array([dN7_dΞ, dN7_dζ])
    dN8_dxy = Jinv @ array([dN8_dΞ, dN8_dζ])

    # ε = B ue
    B = zeros((3, 18))
    B[0, 0] = dN0_dxy[0]
    B[1, 1] = dN0_dxy[1]
    B[2, 0] = dN0_dxy[1]
    B[2, 1] = dN0_dxy[0]
    B[0, 2] = dN1_dxy[0]
    B[1, 3] = dN1_dxy[1]
    B[2, 2] = dN1_dxy[1]
    B[2, 3] = dN1_dxy[0]
    B[0, 4] = dN2_dxy[0]
    B[1, 5] = dN2_dxy[1]
    B[2, 4] = dN2_dxy[1]
    B[2, 5] = dN2_dxy[0]
    B[0, 6] = dN3_dxy[0]
    B[1, 7] = dN3_dxy[1]
    B[2, 6] = dN3_dxy[1]
    B[2, 7] = dN3_dxy[0]
    B[0, 8] = dN4_dxy[0]
    B[1, 9] = dN4_dxy[1]
    B[2, 8] = dN4_dxy[1]
    B[2, 9] = dN4_dxy[0]
    B[0, 10] = dN5_dxy[0]
    B[1, 11] = dN5_dxy[1]
    B[2, 10] = dN5_dxy[1]
    B[2, 11] = dN5_dxy[0]
    B[0, 12] = dN6_dxy[0]
    B[1, 13] = dN6_dxy[1]
    B[2, 12] = dN6_dxy[1]
    B[2, 13] = dN6_dxy[0]
    B[0, 14] = dN7_dxy[0]
    B[1, 15] = dN7_dxy[1]
    B[2, 14] = dN7_dxy[1]
    B[2, 15] = dN7_dxy[0]
    B[0, 16] = dN8_dxy[0]
    B[1, 17] = dN8_dxy[1]
    B[2, 16] = dN8_dxy[1]
    B[2, 17] = dN8_dxy[0]

    u_e = lil_matrix(u_e)
    B = lil_matrix(B)

    #print(B)
    #print(u_e)

    B = B.todense()
    u_e = u_e.toarray()

    #print(B)
    #print(u_e)

    εn = B @ u_e.T
    σn = Eσ @ εn

    return εn, σn

# xy = array([
# [-1,-1],
# [1,-1],
# [1,1],
# [-1,1],
# 	])

# xy = array([
# [0,0],
# [1,0],
# [1,1],
# [0,1],
# 	])

# properties = {}
# properties["E"] = 1.
# properties["nu"] = 0.25
# properties["bx"] = 0
# properties["by"] = 1.
# properties["t"] = 1.

# ke, fe = quad4(xy, properties)

# print(f"ke = {ke}")

# fixed_dofs = [0, 1, 2, 3]
# free_dofs = [4, 5, 6, 7]

# ke_ff = ke[ix_(free_dofs, free_dofs)]
# fe_ff = array([0, -1, 0, -1])

# print(f"ke_ff = {ke_ff}")

# from scipy.linalg import solve

# u = zeros((8,1))
# uf = solve(ke_ff, fe_ff)

# print(f"uf = {uf}")