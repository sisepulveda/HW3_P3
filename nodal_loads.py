import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve

def nodal_loads(xy, properties_load):
    tx = properties_load["tx"]
    e = properties_load["t"]
    ty = properties_load["ty"]

    distx = -(xy[0,0]-xy[1,0])
    disty = -(xy[0,1]-xy[1,1])
    print(f"e = {e} ; disty = {disty} ; ty = {ty}")
    print(f"e = {e} ; distx = {distx} ; tx = {tx}")

    c1x = (disty/2)*e*tx
    c1y = (distx/2)*e*ty
    c2x = (disty/2)*e*tx
    c2y = (distx/2)*e*ty

    fe = [c1x, c1y, c2x, c2y]
    print(fe)
    return fe

def nodal_loads9(xy, properties_load):
    tx = properties_load["tx"]
    e = properties_load["t"]
    ty = properties_load["ty"]

    distx = -(xy[0,0]-xy[1,0])
    disty = -(xy[0,1]-xy[1,1])
    print(f"e = {e} ; disty = {disty} ; ty = {ty}")
    print(f"e = {e} ; distx = {distx} ; tx = {tx}")

    c1x = (disty/2)*e*tx
    c1y = (distx/2)*e*ty
    c2x = (disty/2)*e*tx
    c2y = (distx/2)*e*ty

    fe = [c1x, c1y, c2x, c2y]
    print(fe)
    return fe
