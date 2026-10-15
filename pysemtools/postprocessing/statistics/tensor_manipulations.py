#####################################################################################
## Rotation related manipulations
#####################################################################################

#%% 
# Notation based on Introduction to Continuum Mechanics by Lai, Rubin, Krempl, 4th ed., 2010, Elsevier.
# equation 2.16.5 in the book
def define_rotation_tensor_from_vectors(v1,v2,v3):
    import numpy as np
    return np.stack((v1,v2,v3),axis=4)

#%% rotation of 1D arrray
# Notation based on Introduction to Continuum Mechanics by Lai, Rubin, Krempl, 4th ed., 2010, Elsevier.
# equation 2.17.6 in the book
def rotate_1D_tensor(T1,Qtrans,Q):
    import numpy as np
    return np.matmul(Qtrans,T1)

#%%
# Notation based on Introduction to Continuum Mechanics by Lai, Rubin, Krempl, 4th ed., 2010, Elsevier.
# equation 2.18.5 in the book
def rotate_2D_tensor(T2,Qtrans,Q):
    import numpy as np
    return np.matmul( np.matmul(Qtrans,T2) , Q)

#%%
def convert_2Dtensor_from_6c_to_3x3(T2_6c):
    import numpy as np

    T2_9c = T2_6c[ ..., [0,3,4 , 3,1,5 , 4,5,2] ]
    T2_9c = T2_9c.reshape( T2_9c.shape[0:-1]+(3,3) , order="F" )
    return T2_9c

#%%
def convert_1Dtensor_from_3c_to_3x1(T1):
    import numpy as np

    T1_3x1 = T1[..., np.newaxis]
    return T1_3x1

#%%
def convert_2Dtensor_from_3x3_to_6c(T2_3x3):
    import numpy as np

    T2_6c = T2_3x3.reshape( T2_3x3.shape[:-2]+(9,) , order="F" )
    T2_6c = T2_6c[...,[0,4,8,1,2,5]]
    return T2_6c

#%% helper functions defined for convenience
def rotate_6c2D_tensor(T2,Qtrans,Q):
    T2_rot = convert_2Dtensor_from_6c_to_3x3(T2)
    T2_rot = rotate_2D_tensor(T2_rot,Qtrans,Q)
    T2_rot = convert_2Dtensor_from_3x3_to_6c(T2_rot)
    return T2_rot

def rotate_vector(V,Qtrans,Q):
    V_rot = convert_1Dtensor_from_3c_to_3x1(V)
    V_rot = rotate_1D_tensor(V_rot,Qtrans,Q)
    V_rot = V_rot[...,0]
    return V_rot

#%% cylindrical basis about the z-axis
def define_cylindrical_basis_vectors(XYZ):
    """
    Point-wise orthonormal cylindrical basis (r, theta, z) about the z-axis.

    XYZ has shape (..., 3). Returns v1 = r_hat, v2 = theta_hat = z_hat x r_hat, v3 = z_hat,
    each of shape (..., 3), i.e. rotated components are ordered [U_r, U_theta, U_z].
    """
    import numpy as np

    v1 = np.array(XYZ, dtype=np.float64, copy=True)
    v1[..., 2] = 0.0
    v1 = v1 / np.sqrt(np.sum(v1**2, axis=-1, keepdims=True))

    v3 = np.zeros_like(v1)
    v3[..., 2] = 1.0

    v2 = np.cross(v3, v1, axis=-1)

    return v1, v2, v3

#%% rotation of a full tensor of arbitrary order
# T'_{ij...} = Q_{ai} Q_{bj} ... T_{ab...}, done one index at a time to limit memory use
def rotate_full_tensor(T, Q, order):
    import numpy as np

    letters = "abcdefgh"[:order]
    for p in range(order):
        src = letters[:p] + "z" + letters[p+1:]
        T = np.einsum(f"...z{letters[p]},...{src}->...{letters}", Q, T)
    return T

#%% 3rd order fully symmetric tensors stored with 10 components
# order: UUU, VVV, WWW, UUV, UUW, UVV, VVW, UWW, VWW, UVW
# (i.e. the order of UiUjUk after the reordering in read_interpolated_stat_hdf5_fields)
_IND_10C_TO_27 = [0, 3, 4,  3, 5, 9,  4, 9, 7,
                  3, 5, 9,  5, 1, 6,  9, 6, 8,
                  4, 9, 7,  9, 6, 8,  7, 8, 2]
_IND_27_TO_10C = [0, 13, 26, 1, 2, 4, 14, 8, 17, 5]

def convert_3Dtensor_from_10c_to_3x3x3(T3_10c):
    T3 = T3_10c[..., _IND_10C_TO_27]
    return T3.reshape(T3.shape[:-1] + (3, 3, 3))

def convert_3Dtensor_from_3x3x3_to_10c(T3):
    T3 = T3.reshape(T3.shape[:-3] + (27,))
    return T3[..., _IND_27_TO_10C]

def rotate_10c3D_tensor(T3_10c, Qtrans, Q):
    T3_rot = convert_3Dtensor_from_10c_to_3x3x3(T3_10c)
    T3_rot = rotate_full_tensor(T3_rot, Q, 3)
    return convert_3Dtensor_from_3x3x3_to_10c(T3_rot)

#%% non-symmetric 2nd order tensors stored with 9 components, index 3*i+j (e.g. dUi/dxj)
def rotate_9c2D_tensor(T2_9c, Qtrans, Q):
    T2_rot = T2_9c.reshape(T2_9c.shape[:-1] + (3, 3))
    T2_rot = rotate_2D_tensor(T2_rot, Qtrans, Q)
    return T2_rot.reshape(T2_rot.shape[:-2] + (9,))

#%% rotation of the raw (not yet mean-subtracted) statistics stored as *_struct fields
# Each kind describes how the trailing component axis is laid out:
#   "coords"     : coordinates, not rotated
#   "scalar"     : any number of rotation-invariant components
#   "vec"        : 3 components
#   "sym2"       : symmetric 2nd order, 6 components (11,22,33,12,13,23)
#   "gen2"       : 2nd order, 9 components, index 3*i+j
#   "sym3"       : symmetric 3rd order, 10 components (see rotate_10c3D_tensor)
#   "sym2_grad"  : gradient of a "sym2", 18 components, index 3*p+k (p: sym2 component, k: derivative)
#   "sym3_grad"  : gradient of a "sym3", 30 components, index 3*t+k (t: sym3 component, k: derivative)
#   "lap_scalar" : d2/dxk2 of a scalar, 3 components (k)
#   "lap_vec"    : d2/dxk2 of a vector, 9 components, index 3*i+k
#   "lap_sym2"   : d2/dxk2 of a "sym2", 18 components, index 3*p+k
#   "sdiss"      : [<dS/dxk dS/dxk>, <dS/dxk dUi/dxk>] -> 1 scalar + 3 vector components
#   "ui4"        : diagonal 4th order moments <Ui^4>, cannot be rotated -> NaN
#
# The individual d2/dxk2 are diagonal entries of a Hessian and are not tensor components; only
# their sum over k (the Laplacian) is. The "lap_*" kinds therefore return the rotated Laplacian
# in the k=0 slot and zeros in the k=1,2 slots. Sums over k, which is how the budget routines
# use these fields, are thus preserved; individual slots are not meaningful.
def rotate_stat_struct_field(data, kind, Q):
    import numpy as np

    Qtrans = np.swapaxes(Q, -1, -2)
    lead = data.shape[:-1]

    if kind in ("coords", "scalar"):
        return data

    if kind == "vec":
        return rotate_vector(data, Qtrans, Q)

    if kind == "sym2":
        return rotate_6c2D_tensor(data, Qtrans, Q)

    if kind == "gen2":
        return rotate_9c2D_tensor(data, Qtrans, Q)

    if kind == "sym3":
        return rotate_10c3D_tensor(data, Qtrans, Q)

    if kind == "sym2_grad":
        T = data.reshape(lead + (6, 3))
        T = T[..., [0, 3, 4, 3, 1, 5, 4, 5, 2], :].reshape(lead + (3, 3, 3))
        T = rotate_full_tensor(T, Q, 3)
        T = T.reshape(lead + (9, 3))[..., [0, 4, 8, 1, 2, 5], :]
        return T.reshape(lead + (18,))

    if kind == "sym3_grad":
        T = data.reshape(lead + (10, 3))
        T = T[..., _IND_10C_TO_27, :].reshape(lead + (3, 3, 3, 3))
        T = rotate_full_tensor(T, Q, 4)
        T = T.reshape(lead + (27, 3))[..., _IND_27_TO_10C, :]
        return T.reshape(lead + (30,))

    if kind == "lap_scalar":
        out = np.zeros_like(data)
        out[..., 0] = np.sum(data, axis=-1)
        return out

    if kind == "lap_vec":
        lap = np.sum(data.reshape(lead + (3, 3)), axis=-1)
        out = np.zeros(lead + (3, 3), dtype=data.dtype)
        out[..., 0] = rotate_vector(lap, Qtrans, Q)
        return out.reshape(lead + (9,))

    if kind == "lap_sym2":
        lap = np.sum(data.reshape(lead + (6, 3)), axis=-1)
        out = np.zeros(lead + (6, 3), dtype=data.dtype)
        out[..., 0] = rotate_6c2D_tensor(lap, Qtrans, Q)
        return out.reshape(lead + (18,))

    if kind == "sdiss":
        out = data.copy()
        out[..., 1:4] = rotate_vector(data[..., 1:4], Qtrans, Q)
        return out

    if kind == "ui4":
        return np.full_like(data, np.nan)

    raise ValueError(f"Unknown field kind '{kind}'")
