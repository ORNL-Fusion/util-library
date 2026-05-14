from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from calc_B.circular_coils import CircularCoilSet, build_circular_coil_jackson

try:
    from mpex_modeling_data.coils import (
        get_mpex_currents_by_config as _get_mpex_currents_by_config,
        get_mpex_geometry_table as _get_mpex_geometry_table,
    )
except ImportError as exc:
    _MPEX_MODELING_DATA_IMPORT_ERROR = exc
else:
    _MPEX_MODELING_DATA_IMPORT_ERROR = None


@dataclass(frozen=True)
class MPEXCoilGeometry:
    coil_index: np.ndarray
    coil_name: tuple[str, ...]
    z0: np.ndarray
    cl: np.ndarray
    rr1: np.ndarray
    rr2: np.ndarray
    nturns: np.ndarray
    nlayers: np.ndarray

    @property
    def ncoils(self) -> int:
        return int(self.coil_index.size)

    @property
    def nwind(self) -> np.ndarray:
        return self.nturns * self.nlayers

    @property
    def area(self) -> np.ndarray:
        return self.cl * (self.rr2 - self.rr1)


def get_mpex_geometry_table() -> tuple[list[str], np.ndarray]:
    if _MPEX_MODELING_DATA_IMPORT_ERROR is not None:
        raise ImportError(
            "MPEX coil geometry is provided by the private mpex-modeling-data repo. "
            "Install it or add its python directory to PYTHONPATH."
        ) from _MPEX_MODELING_DATA_IMPORT_ERROR

    return _get_mpex_geometry_table()


def define_mpex_coils() -> MPEXCoilGeometry:
    _, table = get_mpex_geometry_table()
    return MPEXCoilGeometry(
        coil_index=np.asarray(table[:, 0], dtype=int),
        coil_name=tuple(str(name).strip() for name in table[:, 1]),
        z0=np.asarray(table[:, 3], dtype=float),
        cl=np.asarray(table[:, 4], dtype=float),
        rr1=np.asarray(table[:, 5], dtype=float),
        rr2=np.asarray(table[:, 6], dtype=float),
        nturns=np.asarray(table[:, 7], dtype=int),
        nlayers=np.asarray(table[:, 8], dtype=int),
    )


def get_mpex_currents_by_config() -> tuple[tuple[str, ...], np.ndarray]:
    if _MPEX_MODELING_DATA_IMPORT_ERROR is not None:
        raise ImportError(
            "MPEX current configurations are provided by the private mpex-modeling-data repo. "
            "Install it or add its python directory to PYTHONPATH."
        ) from _MPEX_MODELING_DATA_IMPORT_ERROR

    return _get_mpex_currents_by_config()


def setup_mpex_coils(
    current_in: str | np.ndarray,
    verbose: int = 0,
) -> tuple[MPEXCoilGeometry, np.ndarray]:
    coil_geometry = define_mpex_coils()
    config_name, data = get_mpex_currents_by_config()

    if isinstance(current_in, str):
        try:
            ind = config_name.index(current_in)
        except ValueError as exc:
            raise ValueError(f'Unknown MPEX current configuration "{current_in}"') from exc
        current_per_winding = data[:, ind].copy()
    else:
        current_per_winding = np.asarray(current_in, dtype=float).reshape(-1)
        if current_per_winding.size != coil_geometry.ncoils:
            raise ValueError(
                "current_in must be a config name or a numeric array of length "
                f"{coil_geometry.ncoils}"
            )

    if verbose:
        print("-" * 96)
        print("Coil currents set to:")
        print(" ".join(f"{i:8d}" for i in coil_geometry.coil_index))
        print(" ".join(f"{cur:8.1f}" for cur in current_per_winding))
        if isinstance(current_in, str):
            print(f"Using configuration {current_in}")
        print("-" * 96)

    return coil_geometry, current_per_winding


def build_mpex_coils_jackson(
    current_in: str | np.ndarray,
    verbose: int = 0,
) -> tuple[CircularCoilSet, np.ndarray, MPEXCoilGeometry, np.ndarray]:
    coil_geometry, current_per_winding = setup_mpex_coils(current_in, verbose=verbose)

    nwind_tot = int(np.sum(coil_geometry.nwind))
    winding_current = np.zeros(nwind_tot, dtype=float)
    rwind = np.zeros(nwind_tot, dtype=float)
    zwind = np.zeros(nwind_tot, dtype=float)

    i0 = 0
    for i in range(coil_geometry.ncoils):
        coil_an = build_circular_coil_jackson(
            coil_geometry.rr1[i],
            coil_geometry.rr2[i],
            coil_geometry.z0[i],
            coil_geometry.cl[i],
            int(coil_geometry.nturns[i]),
            int(coil_geometry.nlayers[i]),
        )
        nwind_i = int(coil_geometry.nwind[i])
        i1 = i0 + nwind_i
        rwind[i0:i1] = coil_an.rwind
        zwind[i0:i1] = coil_an.zwind
        winding_current[i0:i1] = current_per_winding[i]
        i0 = i1

    return (
        CircularCoilSet(rwind=rwind, zwind=zwind),
        winding_current,
        coil_geometry,
        current_per_winding,
    )


def build_mpex_coils_jackson_hybrid(
    current_in: str | np.ndarray,
    simplify_coils: None | np.ndarray | list[int] | list[bool] = None,
    nturns_simple: int = 3,
    nlayers_simple: int = 2,
    verbose: int = 0,
) -> tuple[CircularCoilSet, np.ndarray, MPEXCoilGeometry, np.ndarray]:
    if nturns_simple <= 0 or nlayers_simple <= 0:
        raise ValueError("nturns_simple and nlayers_simple must be positive")
    if int(nturns_simple) != nturns_simple or int(nlayers_simple) != nlayers_simple:
        raise ValueError("nturns_simple and nlayers_simple must be integers")

    coil_geometry, current_per_winding = setup_mpex_coils(current_in, verbose=verbose)
    simplify_mask = _get_simplify_mask(coil_geometry, simplify_coils)
    nsimple = int(nturns_simple * nlayers_simple)

    nwind_tot = int(np.sum(coil_geometry.nwind[~simplify_mask]) + nsimple * np.sum(simplify_mask))
    winding_current = np.zeros(nwind_tot, dtype=float)
    rwind = np.zeros(nwind_tot, dtype=float)
    zwind = np.zeros(nwind_tot, dtype=float)

    i0 = 0
    for i in range(coil_geometry.ncoils):
        if simplify_mask[i]:
            coil_an = build_circular_coil_jackson(
                coil_geometry.rr1[i],
                coil_geometry.rr2[i],
                coil_geometry.z0[i],
                coil_geometry.cl[i],
                int(nturns_simple),
                int(nlayers_simple),
            )
            nwind_i = nsimple
            cur_i = current_per_winding[i] * coil_geometry.nwind[i] / nsimple
        else:
            coil_an = build_circular_coil_jackson(
                coil_geometry.rr1[i],
                coil_geometry.rr2[i],
                coil_geometry.z0[i],
                coil_geometry.cl[i],
                int(coil_geometry.nturns[i]),
                int(coil_geometry.nlayers[i]),
            )
            nwind_i = int(coil_geometry.nwind[i])
            cur_i = current_per_winding[i]

        i1 = i0 + nwind_i
        rwind[i0:i1] = coil_an.rwind
        zwind[i0:i1] = coil_an.zwind
        winding_current[i0:i1] = cur_i
        i0 = i1

    if verbose:
        print(
            f"Hybrid builder simplified {int(np.sum(simplify_mask))} "
            f"of {coil_geometry.ncoils} coils."
        )
        print(f"Simplified coils use {nturns_simple} turns and {nlayers_simple} layers.")
        if np.any(simplify_mask):
            simplified = np.asarray(coil_geometry.coil_index)[simplify_mask]
            print("Simplified coil indices:", " ".join(str(int(v)) for v in simplified))

    return (
        CircularCoilSet(rwind=rwind, zwind=zwind),
        winding_current,
        coil_geometry,
        current_per_winding,
    )


def get_coil_cross_sections(coil_geometry: MPEXCoilGeometry) -> tuple[np.ndarray, np.ndarray]:
    rcoil = np.column_stack(
        [coil_geometry.rr1, coil_geometry.rr2, coil_geometry.rr2, coil_geometry.rr1, coil_geometry.rr1]
    )
    zcoil = np.column_stack(
        [
            coil_geometry.z0,
            coil_geometry.z0,
            coil_geometry.z0 + coil_geometry.cl,
            coil_geometry.z0 + coil_geometry.cl,
            coil_geometry.z0,
        ]
    )
    return rcoil, zcoil


def get_mpex_plot_defaults() -> dict[str, object]:
    return {
        "config_name": ["D3-6", "D2-2", "D1-1"],
        "field_models": ["jackson", "hybrid"],
        "primary_model": "hybrid",
        "hybrid_simplify_coils": None,
        "hybrid_nturns": 3,
        "hybrid_nlayers": 3,
        "n_test": 200,
        "z_limits": (-2.0, 10.0),
        "plot_coils": False,
        "verbose": 1,
    }


def _get_simplify_mask(
    coil_geometry: MPEXCoilGeometry,
    simplify_coils: None | np.ndarray | list[int] | list[bool],
) -> np.ndarray:
    ncoils = coil_geometry.ncoils

    if simplify_coils is None:
        rmid = 0.5 * (coil_geometry.rr1 + coil_geometry.rr2)
        return rmid > 0.4

    simplify_arr = np.asarray(simplify_coils)
    if simplify_arr.dtype == bool:
        if simplify_arr.size != ncoils:
            raise ValueError(f"Logical simplify_coils array must have length {ncoils}")
        return simplify_arr.reshape(ncoils).astype(bool)

    simplify_idx = simplify_arr.astype(int).reshape(-1)
    if np.any(simplify_idx < 1) or np.any(simplify_idx > ncoils):
        raise ValueError(f"simplify_coils indices must be between 1 and {ncoils}")

    simplify_mask = np.zeros(ncoils, dtype=bool)
    simplify_mask[simplify_idx - 1] = True
    return simplify_mask
