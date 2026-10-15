"""
Tests for rotating and averaging the raw interpolated statistics BEFORE the fluctuating quantities
and budgets are formed (rotate_and_average_interpolated_*stat_fields), compared to the legacy
workflow where the budgets are formed point-wise in Cartesian coordinates and only then rotated
and averaged.
"""

import h5py
import numpy as np
import pytest

from pysemtools.postprocessing.statistics.tensor_manipulations import (
    define_cylindrical_basis_vectors,
    define_rotation_tensor_from_vectors,
    rotate_full_tensor,
    rotate_10c3D_tensor,
    rotate_9c2D_tensor,
    rotate_2D_tensor,
    rotate_6c2D_tensor,
    rotate_vector,
    rotate_stat_struct_field,
    convert_3Dtensor_from_10c_to_3x3x3,
    convert_3Dtensor_from_3x3x3_to_10c,
)
from pysemtools.postprocessing.statistics.RS_budgets_interpolatedPoints_notMPI import (
    FLUID_STRUCT_FIELD_KINDS,
    calculate_budgets_in_Cartesian,
    rotate_and_average_interpolated_stat_fields,
)
from pysemtools.postprocessing.statistics.passive_scalar_budgets_interpolatedPoints_notMPI import (
    SCALAR_STRUCT_FIELD_KINDS,
    calculate_scalar_budgets_in_Cartesian,
    rotate_and_average_interpolated_scalar_stat_fields,
)

NZ, NTHETA = 5, 12
R0, L_THETA = 2.0, np.pi / 3

FLUID_NCOMP = {
    "XYZ": 3, "UVW": 3, "UiUj": 6, "UiUjUk": 10, "Ui4": 3, "P": 1, "P2": 1, "P3": 1, "P4": 1,
    "PUi": 3, "PGij": 9, "pseudoDiss": 6, "dUidxj": 9, "d2Uidx2": 9, "dPdx": 3, "d2Pdx2": 3,
    "dPUidxj": 9, "dUiUjdx": 18, "d2UiUjdx2": 18, "dUiUjUkdx": 30,
}
SCALAR_NCOMP = {
    "XYZ": 3, "S": 1, "S2": 1, "S3": 1, "S4": 1, "PS": 1, "UiS": 3, "UiS2": 3, "UiUjS": 6,
    "PdSdxi": 3, "UidSdxj": 9, "SdUidxj": 9, "SDiss": 4, "dSdxj": 3, "d2Sdxj2": 3, "dS2dxj": 3,
    "d2S2dxj2": 3, "dUiSdxj": 9, "dUjS2dxj": 1, "dUiUjSdxj": 3, "dPSdxj": 3, "dUidSdxjdxj": 3,
    "dSdUidxjdxj": 3,
}

# tensor kind of each output of the budget routines (as rotated by the legacy scripts)
FLUID_OUTPUT_KINDS = {
    "UVW": "vec", "Rij": "sym2", "P": "scalar", "P2": "scalar", "P3": "scalar", "P4": "scalar",
    "UiUjUk": "sym3", "dUidXj": "dUidXj", "Momentum_convection": "vec", "Prod_ij": "sym2",
    "TKE_prod": "scalar", "Conv_ij": "sym2", "TKE_conv": "scalar", "Diss_ij": "sym2",
    "TKE_diss": "scalar", "PRS_ij": "sym2", "dpdx": "vec", "PressTrans_ij": "sym2",
    "VelPressGrad_ij": "sym2", "TKE_VelPressGrad": "scalar", "ViscDiff_ij": "sym2",
    "TKE_ViscDiff": "scalar", "Momentum_viscous_diffusion": "vec", "Momentum_pressure": "vec",
    "Momentum_turb_diffusion": "vec", "TurbTrans_ij": "sym2", "TKE_turbTrans": "scalar",
    "Momentum_residual": "vec", "TKE_residual": "scalar",
}
# not tensorial, so not comparable: Ui4 (only diagonal 4th moments), Momentum_convectionRes
FLUID_SKIPPED_OUTPUTS = {"Ui4", "Momentum_convectionRes", "XYZ", "Rer_here"}

SCALAR_OUTPUT_KINDS = {
    "S": "scalar", "S2": "scalar", "S3": "scalar", "S4": "scalar", "UiS": "vec", "UiS2": "vec",
    "UiUjS": "sym2", "PS": "scalar", "dSdXj": "vec",
    "S_convection": "scalar", "S_turb_diffusion": "scalar", "S_molecular_diffusion": "scalar",
    "S_residual": "scalar", "S2_convection": "scalar", "S2_production": "scalar",
    "S2_turb_diffusion": "scalar", "S2_molecular_diffusion": "scalar", "S2_dissipation": "scalar",
    "S2_residual": "scalar", "UiS_convection": "vec", "UiS_velocity_gradient_production": "vec",
    "UiS_scalar_gradient_production": "vec", "UiS_turb_diffusion": "vec",
    "UiS_molecular_diffusion": "vec", "UiS_scalar_pressure_correlation_gradient": "vec",
    "UiS_pressure_scalar_gradient_correlation": "vec", "UiS_dissipation": "vec",
    "UiS_residual": "vec",
}
SCALAR_SKIPPED_OUTPUTS = {"XYZ", "Rer_here", "Prt_here"}


def theta_mean(x):
    return np.mean(x, axis=2, keepdims=True)


def identity(x):
    return x


###########################################################################################
# legacy rotation (copied from the original disk extraction scripts)
###########################################################################################
def legacy_rotate_10c3D_tensor(T3_10c, Qtrans, Q):
    indices_map = [0, 3, 4, 3, 5, 9, 4, 9, 7, 3, 5, 9, 5, 1, 6, 9, 6, 8, 4, 9, 7, 9, 6, 8, 7, 8, 2]
    T3_full = T3_10c[..., indices_map]
    T3_full = T3_full.reshape(T3_full.shape[:-1] + (3, 3, 3))
    T3_rot = np.einsum("...im,...jn,...kp,...mnp->...ijk", Qtrans, Qtrans, Qtrans, T3_full)
    return np.stack(
        [T3_rot[..., 0, 0, 0], T3_rot[..., 1, 1, 1], T3_rot[..., 2, 2, 2],
         T3_rot[..., 0, 0, 1], T3_rot[..., 0, 0, 2], T3_rot[..., 0, 1, 1],
         T3_rot[..., 1, 1, 2], T3_rot[..., 0, 2, 2], T3_rot[..., 1, 2, 2],
         T3_rot[..., 0, 1, 2]], axis=-1)


def legacy_rotate_output(data, kind, Q):
    Qtrans = np.swapaxes(Q, -1, -2)
    if kind == "scalar":
        return data
    if kind == "vec":
        return rotate_vector(data, Qtrans, Q)
    if kind == "sym2":
        return rotate_6c2D_tensor(data, Qtrans, Q)
    if kind == "sym3":
        return legacy_rotate_10c3D_tensor(data, Qtrans, Q)
    if kind == "dUidXj":
        # as in the scripts: Fortran-order reshape, rotate, keep as 3x3
        T = np.reshape(data, data.shape[:-1] + (3, 3), order="F")
        return rotate_2D_tensor(T, Qtrans, Q)
    raise ValueError(kind)


def new_output(data, kind):
    if kind == "dUidXj":
        return np.reshape(data, data.shape[:-1] + (3, 3), order="F")
    return data


###########################################################################################
# data generation
###########################################################################################
def make_xyz():
    theta_full = np.linspace(0.0, L_THETA, NTHETA + 1)
    theta = 0.5 * (theta_full[1:] + theta_full[:-1])
    z = np.linspace(0.0, 1.0, NZ)
    Z, T = np.meshgrid(z, theta, indexing="ij")
    XYZ = np.stack([R0 * np.cos(T), R0 * np.sin(T), Z], axis=-1)
    return XYZ[np.newaxis, ...]  # (1, NZ, NTHETA, 3)


def cylindrical_Q(XYZ):
    return define_rotation_tensor_from_vectors(*define_cylindrical_basis_vectors(XYZ))


def random_raw_fields(rng, ncomp, XYZ, homogeneous):
    """
    Random raw fields in Cartesian components. If homogeneous, the cylindrical components only
    depend on z (constant along theta), as for azimuthally homogeneous statistics.
    """
    Q = cylindrical_Q(XYZ)
    kinds = {**FLUID_STRUCT_FIELD_KINDS, **SCALAR_STRUCT_FIELD_KINDS}
    fields = {}
    for name, nc in ncomp.items():
        if name == "XYZ":
            fields[name] = XYZ.copy()
            continue
        if homogeneous:
            cyl = rng.uniform(-1, 1, size=(1, NZ, 1, nc)) * np.ones((1, NZ, NTHETA, nc))
            kind = kinds[name]
            if kind == "ui4":
                fields[name] = cyl
            else:
                # cylindrical -> Cartesian is the rotation with Q^T
                fields[name] = rotate_stat_struct_field(cyl, kind, np.swapaxes(Q, -1, -2))
        else:
            fields[name] = rng.uniform(-1, 1, size=(1, NZ, NTHETA, nc))
    return fields


def write_struct(path, fields, extra):
    with h5py.File(path, "w") as hf:
        for name, data in fields.items():
            hf.create_dataset(name + "_struct", data=data)
        for name, val in extra.items():
            hf.create_dataset(name, data=val)


def read_all(path):
    with h5py.File(path, "r") as hf:
        return {k: np.array(hf[k]) for k in hf.keys()}


###########################################################################################
# the two workflows
###########################################################################################
def run_fluid_workflows(tmp_path, fluid_fields, av_func):
    d = str(tmp_path)
    write_struct(tmp_path / "raw_fluid.hdf5", fluid_fields, {"Rer_here": 1000.0})

    # legacy: budgets point-wise in Cartesian -> rotate -> average
    calculate_budgets_in_Cartesian(d, "raw_fluid.hdf5", "legacy_fluid_budgets.hdf5")
    legacy = read_all(tmp_path / "legacy_fluid_budgets.hdf5")
    Q = cylindrical_Q(legacy["XYZ"])
    legacy = {k: av_func(legacy_rotate_output(legacy[k], kind, Q))
              for k, kind in FLUID_OUTPUT_KINDS.items()}

    # new: rotate raw -> average -> budgets
    rotate_and_average_interpolated_stat_fields(
        d, "raw_fluid.hdf5", "cyl_fluid.hdf5",
        basis_func=define_cylindrical_basis_vectors, av_func=av_func)
    calculate_budgets_in_Cartesian(d, "cyl_fluid.hdf5", "new_fluid_budgets.hdf5")
    new = read_all(tmp_path / "new_fluid_budgets.hdf5")
    assert set(new) == set(FLUID_OUTPUT_KINDS) | FLUID_SKIPPED_OUTPUTS
    new = {k: new_output(new[k], kind) for k, kind in FLUID_OUTPUT_KINDS.items()}
    return legacy, new


def run_scalar_workflows(tmp_path, fluid_fields, scalar_fields, av_func):
    d = str(tmp_path)
    write_struct(tmp_path / "raw_fluid.hdf5", fluid_fields, {"Rer_here": 1000.0})
    write_struct(tmp_path / "raw_scalar.hdf5", scalar_fields,
                 {"Rer_here": 1000.0, "Prt_here": 0.71})

    calculate_scalar_budgets_in_Cartesian(
        d, "raw_scalar.hdf5", "raw_fluid.hdf5", "legacy_scalar_budgets.hdf5")
    legacy = read_all(tmp_path / "legacy_scalar_budgets.hdf5")
    Q = cylindrical_Q(legacy["XYZ"])
    legacy = {k: av_func(legacy_rotate_output(legacy[k], kind, Q))
              for k, kind in SCALAR_OUTPUT_KINDS.items()}

    rotate_and_average_interpolated_stat_fields(
        d, "raw_fluid.hdf5", "cyl_fluid.hdf5",
        basis_func=define_cylindrical_basis_vectors, av_func=av_func)
    rotate_and_average_interpolated_scalar_stat_fields(
        d, "raw_scalar.hdf5", "cyl_scalar.hdf5",
        basis_func=define_cylindrical_basis_vectors, av_func=av_func)
    calculate_scalar_budgets_in_Cartesian(
        d, "cyl_scalar.hdf5", "cyl_fluid.hdf5", "new_scalar_budgets.hdf5")
    new = read_all(tmp_path / "new_scalar_budgets.hdf5")
    assert set(new) == set(SCALAR_OUTPUT_KINDS) | SCALAR_SKIPPED_OUTPUTS
    new = {k: new[k] for k in SCALAR_OUTPUT_KINDS}
    return legacy, new


def assert_all_close(legacy, new, rtol=1e-11):
    for k in legacy:
        assert legacy[k].shape == new[k].shape, k
        scale = max(np.max(np.abs(legacy[k])), 1.0)
        np.testing.assert_allclose(new[k], legacy[k], rtol=0, atol=rtol * scale, err_msg=k)


###########################################################################################
# tests of the rotation helpers
###########################################################################################
def test_10c_roundtrip_and_symmetry():
    rng = np.random.default_rng(0)
    T10 = rng.normal(size=(4, 10))
    T27 = convert_3Dtensor_from_10c_to_3x3x3(T10)
    for perm in [(0, 2, 1), (1, 0, 2), (2, 1, 0), (1, 2, 0), (2, 0, 1)]:
        np.testing.assert_array_equal(T27, np.transpose(T27, (0,) + tuple(p + 1 for p in perm)))
    np.testing.assert_array_equal(convert_3Dtensor_from_3x3x3_to_10c(T27), T10)


def test_rotation_helpers_against_references():
    rng = np.random.default_rng(1)
    Q = cylindrical_Q(make_xyz())
    Qtrans = np.swapaxes(Q, -1, -2)
    shape = Q.shape[:-2]

    T10 = rng.normal(size=shape + (10,))
    np.testing.assert_allclose(rotate_10c3D_tensor(T10, Qtrans, Q),
                               legacy_rotate_10c3D_tensor(T10, Qtrans, Q), atol=1e-14)

    T9 = rng.normal(size=shape + (9,))
    ref = rotate_2D_tensor(T9.reshape(shape + (3, 3)), Qtrans, Q).reshape(shape + (9,))
    np.testing.assert_allclose(rotate_9c2D_tensor(T9, Qtrans, Q), ref, atol=1e-14)
    # the Fortran-order reshape used in the scripts gives the transpose, which rotates the same
    refF = rotate_2D_tensor(T9.reshape(shape + (3, 3), order="F"), Qtrans, Q)
    np.testing.assert_allclose(
        rotate_9c2D_tensor(T9, Qtrans, Q).reshape(shape + (3, 3), order="F"), refF, atol=1e-14)

    T4 = rng.normal(size=shape + (3, 3, 3, 3))
    back = rotate_full_tensor(rotate_full_tensor(T4, Q, 4), Qtrans, 4)
    np.testing.assert_allclose(back, T4, atol=1e-13)
    ref4 = np.einsum("...ai,...bj,...ck,...dl,...abcd->...ijkl", Q, Q, Q, Q, T4)
    np.testing.assert_allclose(rotate_full_tensor(T4, Q, 4), ref4, atol=1e-13)


def test_gradient_kinds_match_full_tensor_rotation():
    rng = np.random.default_rng(2)
    Q = cylindrical_Q(make_xyz())
    shape = Q.shape[:-2]

    # sym2_grad: d<uiuj>/dxk from a symmetric 3x3x3 tensor
    T = rng.normal(size=shape + (3, 3, 3))
    T = T + np.swapaxes(T, -3, -2)
    data = T.reshape(shape + (9, 3))[..., [0, 4, 8, 1, 2, 5], :].reshape(shape + (18,))
    ref = rotate_full_tensor(T, Q, 3).reshape(shape + (9, 3))[..., [0, 4, 8, 1, 2, 5], :]
    np.testing.assert_allclose(rotate_stat_struct_field(data, "sym2_grad", Q),
                               ref.reshape(shape + (18,)), atol=1e-13)

    # sym3_grad: d<uiujuk>/dxl
    T10 = rng.normal(size=shape + (10, 3))
    T = convert_3Dtensor_from_10c_to_3x3x3(np.moveaxis(T10, -1, -2))  # (..., 3(l), 3, 3, 3)
    T = np.moveaxis(T, -4, -1)  # (..., i, j, k, l)
    Tr = rotate_full_tensor(T, Q, 4)
    ref = convert_3Dtensor_from_3x3x3_to_10c(np.moveaxis(Tr, -1, -4))  # (..., l, 10)
    ref = np.moveaxis(ref, -2, -1).reshape(shape + (30,))
    np.testing.assert_allclose(rotate_stat_struct_field(T10.reshape(shape + (30,)), "sym3_grad", Q),
                               ref, atol=1e-13)

    # Laplacian kinds keep the sum over the derivative direction
    lap = rng.normal(size=shape + (9,))
    out = rotate_stat_struct_field(lap, "lap_vec", Q).reshape(shape + (3, 3))
    np.testing.assert_allclose(out.sum(-1),
                               rotate_vector(lap.reshape(shape + (3, 3)).sum(-1),
                                             np.swapaxes(Q, -1, -2), Q), atol=1e-14)


###########################################################################################
# tests of the workflows
###########################################################################################
@pytest.mark.parametrize("homogeneous", [False, True])
def test_fluid_pointwise_covariance(tmp_path, homogeneous):
    """Without averaging, rotating before or after the budgets must give the same result."""
    rng = np.random.default_rng(10)
    fluid = random_raw_fields(rng, FLUID_NCOMP, make_xyz(), homogeneous)
    legacy, new = run_fluid_workflows(tmp_path, fluid, identity)
    assert_all_close(legacy, new)


def test_fluid_homogeneous_equivalence(tmp_path):
    """Azimuthally homogeneous statistics: the new workflow reproduces the old one."""
    rng = np.random.default_rng(11)
    fluid = random_raw_fields(rng, FLUID_NCOMP, make_xyz(), homogeneous=True)
    legacy, new = run_fluid_workflows(tmp_path, fluid, theta_mean)
    assert_all_close(legacy, new)


def test_fluid_inhomogeneous_dispersive_stress(tmp_path):
    """With theta-dependent means, the difference is the dispersive stress <U~i U~j>_theta."""
    rng = np.random.default_rng(12)
    XYZ = make_xyz()
    fluid = random_raw_fields(rng, FLUID_NCOMP, XYZ, homogeneous=False)
    legacy, new = run_fluid_workflows(tmp_path, fluid, theta_mean)

    Q = cylindrical_Q(XYZ)
    U = rotate_vector(fluid["UVW"], np.swapaxes(Q, -1, -2), Q)
    Ut = U - theta_mean(U)
    i, j = [0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]
    dispersive = theta_mean(Ut[..., i] * Ut[..., j])
    Pt = fluid["P"] - theta_mean(fluid["P"])

    np.testing.assert_allclose(new["UVW"], legacy["UVW"], atol=1e-14)
    np.testing.assert_allclose(new["Rij"] - legacy["Rij"], dispersive, atol=1e-13)
    np.testing.assert_allclose(new["P2"] - legacy["P2"], theta_mean(Pt**2), atol=1e-13)
    assert np.max(np.abs(dispersive)) > 1e-2


@pytest.mark.parametrize("homogeneous", [False, True])
def test_scalar_pointwise_covariance(tmp_path, homogeneous):
    rng = np.random.default_rng(20)
    XYZ = make_xyz()
    fluid = random_raw_fields(rng, FLUID_NCOMP, XYZ, homogeneous)
    scalar = random_raw_fields(rng, SCALAR_NCOMP, XYZ, homogeneous)
    legacy, new = run_scalar_workflows(tmp_path, fluid, scalar, identity)
    assert_all_close(legacy, new)


def test_scalar_homogeneous_equivalence(tmp_path):
    rng = np.random.default_rng(21)
    XYZ = make_xyz()
    fluid = random_raw_fields(rng, FLUID_NCOMP, XYZ, homogeneous=True)
    scalar = random_raw_fields(rng, SCALAR_NCOMP, XYZ, homogeneous=True)
    legacy, new = run_scalar_workflows(tmp_path, fluid, scalar, theta_mean)
    assert_all_close(legacy, new)


def test_scalar_inhomogeneous_dispersive_flux(tmp_path):
    rng = np.random.default_rng(22)
    XYZ = make_xyz()
    fluid = random_raw_fields(rng, FLUID_NCOMP, XYZ, homogeneous=False)
    scalar = random_raw_fields(rng, SCALAR_NCOMP, XYZ, homogeneous=False)
    legacy, new = run_scalar_workflows(tmp_path, fluid, scalar, theta_mean)

    Q = cylindrical_Q(XYZ)
    U = rotate_vector(fluid["UVW"], np.swapaxes(Q, -1, -2), Q)
    S = scalar["S"]
    dispersive = theta_mean((U - theta_mean(U)) * (S - theta_mean(S)))
    np.testing.assert_allclose(new["UiS"] - legacy["UiS"], dispersive, atol=1e-13)
    np.testing.assert_allclose(new["S2"] - legacy["S2"],
                               theta_mean((S - theta_mean(S))**2), atol=1e-13)


def test_unknown_field_raises(tmp_path):
    rng = np.random.default_rng(30)
    fluid = random_raw_fields(rng, FLUID_NCOMP, make_xyz(), homogeneous=False)
    fluid["mystery"] = rng.normal(size=(1, NZ, NTHETA, 3))
    write_struct(tmp_path / "raw.hdf5", fluid, {"Rer_here": 1000.0})
    with pytest.raises(KeyError):
        rotate_and_average_interpolated_stat_fields(
            str(tmp_path), "raw.hdf5", "out.hdf5",
            basis_func=define_cylindrical_basis_vectors, av_func=theta_mean)
