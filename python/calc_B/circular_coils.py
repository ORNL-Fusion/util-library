from __future__ import annotations

from dataclasses import dataclass

import numpy as np


MU0_OVER_4PI = 1.0e-7


@dataclass(frozen=True)
class CircularCoilSet:
    rwind: np.ndarray
    zwind: np.ndarray

    def __post_init__(self) -> None:
        if self.rwind.shape != self.zwind.shape:
            raise ValueError("rwind and zwind must have the same shape")


def build_circular_coil_jackson(
    r1: float,
    r2: float,
    z1: float,
    dz: float,
    nturns: int,
    nlayers: int,
) -> CircularCoilSet:
    """Python equivalent of MATLAB ``build_circular_coil_jackson``."""
    if nturns <= 0 or nlayers <= 0:
        raise ValueError("nturns and nlayers must be positive")

    fw = dz / nturns
    fh = (r2 - r1) / nlayers

    zwind = np.tile(np.linspace(z1 + fw / 2.0, z1 + dz - fw / 2.0, nturns), nlayers)
    rwind = np.tile(np.linspace(r1 + fh / 2.0, r2 - fh / 2.0, nlayers), nturns)
    return CircularCoilSet(rwind=rwind, zwind=zwind)


def bfield_circular_coils_vectorized(
    coil: CircularCoilSet,
    current: np.ndarray,
    r: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Python equivalent of MATLAB ``bfield_circular_coils_vectorized``.

    The on-axis path is fully self-contained. Off-axis evaluation uses
    complete elliptic integrals and therefore requires SciPy.
    """
    r_arr = np.asarray(r, dtype=float)
    z_arr = np.asarray(z, dtype=float)
    cur = np.asarray(current, dtype=float)
    awind = np.asarray(coil.rwind, dtype=float).reshape(1, -1)
    dwind = np.asarray(coil.zwind, dtype=float).reshape(1, -1)

    if r_arr.shape != z_arr.shape:
        raise ValueError("r and z must have the same shape")
    if awind.size != cur.size:
        raise ValueError("coil and current must have the same number of windings")

    npts = r_arr.size
    rvec = r_arr.reshape(-1)
    zvec = z_arr.reshape(-1)
    cur_row = cur.reshape(1, -1)

    br_sum = np.zeros(npts, dtype=float)
    bz_sum = np.zeros(npts, dtype=float)

    is_axis = np.isclose(rvec, 0.0)
    is_off_axis = ~is_axis

    if np.any(is_axis):
        z0 = zvec[is_axis].reshape(-1, 1)
        dz0 = z0 - dwind
        bz_axis = 2.0 * np.pi * MU0_OVER_4PI * awind**2 / (awind**2 + dz0**2) ** 1.5
        bz_sum[is_axis] = (bz_axis * cur_row).sum(axis=1)

    if np.any(is_off_axis):
        try:
            from scipy.special import ellipe, ellipk
        except ImportError as exc:
            raise ImportError(
                "Off-axis circular-coil evaluation requires scipy.special. "
                "On-axis evaluation works without SciPy."
            ) from exc

        ro = rvec[is_off_axis].reshape(-1, 1)
        zo = zvec[is_off_axis].reshape(-1, 1)
        dz = zo - dwind
        m = 4.0 * (ro * awind) / ((ro + awind) ** 2 + dz**2)
        sm = np.sqrt(m)
        k_comp = ellipk(m)
        e_comp = ellipe(m)

        sa = np.sqrt(awind)
        sroa = np.sqrt(ro / awind)
        m1 = m - 1.0

        br_terms = (
            -2.0
            * MU0_OVER_4PI
            / ro**1.5
            * dz
            / sa
            * (sm / 2.0 * k_comp - sm / 4.0 * (m - 2.0) / m1 * e_comp)
        )
        bz_terms = (
            2.0
            * MU0_OVER_4PI
            / ro
            * (
                sroa * sm / 2.0 * k_comp
                - sm / (4.0 * m1) * (sroa * (m - 2.0) + m / sroa) * e_comp
            )
        )

        br_sum[is_off_axis] = (br_terms * cur_row).sum(axis=1)
        bz_sum[is_off_axis] = (bz_terms * cur_row).sum(axis=1)

    return br_sum.reshape(r_arr.shape), bz_sum.reshape(r_arr.shape)
