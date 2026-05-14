from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Allow this file to be run directly from the source tree, while keeping
# the normal import path used by setup_paths.sh and editable installs.
if __package__ in (None, ""):
    python_root = Path(__file__).resolve().parents[2]
    if str(python_root) not in sys.path:
        sys.path.insert(0, str(python_root))

from calc_B.circular_coils import bfield_circular_coils_vectorized
from MPEX.COILS.MPEX.coils import (
    build_mpex_coils_jackson,
    build_mpex_coils_jackson_hybrid,
    get_coil_cross_sections,
    get_mpex_plot_defaults,
)


def plot_bfield_on_axis_mpex(
    config_name: list[str] | tuple[str, ...] | None = None,
    field_models: list[str] | tuple[str, ...] | None = None,
    primary_model: str | None = None,
    hybrid_simplify_coils=None,
    hybrid_nturns: int | None = None,
    hybrid_nlayers: int | None = None,
    n_test: int | None = None,
    z_limits: tuple[float, float] | None = None,
    plot_coils: bool | None = None,
    verbose: int | None = None,
    ax=None,
    coil_ax=None,
):
    """Reproduce the MATLAB ``plot_bfield_on_axis_MPEX.m`` figure in Python."""
    defaults = get_mpex_plot_defaults()
    config_name = list(defaults["config_name"] if config_name is None else config_name)
    field_models = list(defaults["field_models"] if field_models is None else field_models)
    primary_model = defaults["primary_model"] if primary_model is None else primary_model
    hybrid_simplify_coils = (
        defaults["hybrid_simplify_coils"]
        if hybrid_simplify_coils is None
        else hybrid_simplify_coils
    )
    hybrid_nturns = defaults["hybrid_nturns"] if hybrid_nturns is None else hybrid_nturns
    hybrid_nlayers = defaults["hybrid_nlayers"] if hybrid_nlayers is None else hybrid_nlayers
    n_test = defaults["n_test"] if n_test is None else n_test
    z_limits = defaults["z_limits"] if z_limits is None else z_limits
    plot_coils = defaults["plot_coils"] if plot_coils is None else plot_coils
    verbose = defaults["verbose"] if verbose is None else verbose

    if primary_model not in field_models:
        raise ValueError("primary_model must be one of field_models")
    if n_test < 2:
        raise ValueError("n_test must be at least 2")
    if z_limits[0] >= z_limits[1]:
        raise ValueError("z_limits must be ordered as (z_min, z_max)")

    import matplotlib.pyplot as plt

    z_test = np.linspace(z_limits[0], z_limits[1], n_test)
    r_test = np.zeros_like(z_test)

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    else:
        fig = ax.figure

    ax.set_facecolor("white")
    ax.grid(True)
    ax.tick_params(labelsize=12)
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, len(config_name)))

    coil_fig = None
    model_results: dict[str, dict[str, dict[str, np.ndarray]]] = {}

    for i, config in enumerate(config_name):
        model_data: dict[str, dict[str, np.ndarray]] = {}

        for model_name in field_models:
            if model_name == "jackson":
                coil, current, coil_geometry, current_per_winding = build_mpex_coils_jackson(
                    config, verbose=verbose
                )
            elif model_name == "hybrid":
                coil, current, coil_geometry, current_per_winding = build_mpex_coils_jackson_hybrid(
                    config,
                    simplify_coils=hybrid_simplify_coils,
                    nturns_simple=hybrid_nturns,
                    nlayers_simple=hybrid_nlayers,
                    verbose=verbose,
                )
            else:
                raise ValueError(f"Unsupported field model: {model_name}")

            br, bz = bfield_circular_coils_vectorized(coil, current, r_test, z_test)
            model_data[model_name] = {
                "Btot": np.sqrt(br**2 + bz**2),
                "coil_geometry": coil_geometry,
                "current_per_winding": current_per_winding,
            }

        if plot_coils and i == 0:
            rcoil, zcoil = get_coil_cross_sections(model_data[primary_model]["coil_geometry"])
            coil_fig, coil_ax = _plot_coil_cross_section(
                rcoil,
                zcoil,
                current=model_data[primary_model]["current_per_winding"],
                ax=coil_ax,
            )

        for model_name in field_models:
            line_style = "-" if model_name == primary_model else "--"
            ax.plot(
                z_test,
                model_data[model_name]["Btot"],
                color=colors[i],
                linestyle=line_style,
                linewidth=2.0,
                label=f"{config} ({model_name})",
            )

        model_results[config] = model_data

    ax.legend(loc="best")
    ax.set_xlabel("Z (m)", fontsize=13)
    ax.set_ylabel("|B| on-axis (T)", fontsize=13)

    return {
        "fig": fig,
        "ax": ax,
        "coil_fig": coil_fig,
        "coil_ax": coil_ax,
        "z_test": z_test,
        "r_test": r_test,
        "model_results": model_results,
    }


def _plot_coil_cross_section(rcoil, zcoil, current=None, ax=None):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    else:
        fig = ax.figure

    ax.set_facecolor("white")
    ax.grid(True)
    ax.tick_params(labelsize=12)

    cmap = plt.cm.viridis
    norm = None
    if current is not None:
        current = np.asarray(current, dtype=float)
        vmin = float(np.min(current))
        vmax = float(np.max(current))
        norm = plt.Normalize(vmin=vmin, vmax=vmax)

    for i in range(rcoil.shape[0]):
        ax.plot(zcoil[i, :], rcoil[i, :], color="k", linewidth=0.8)
        if current is not None:
            ax.fill(zcoil[i, :], rcoil[i, :], color=cmap(norm(current[i])), alpha=0.85)

    if current is not None:
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label="Current per winding (A)")

    ax.set_xlabel("Z (m)", fontsize=13)
    ax.set_ylabel("R (m)", fontsize=13)
    return fig, ax


def _build_parser() -> argparse.ArgumentParser:
    defaults = get_mpex_plot_defaults()
    parser = argparse.ArgumentParser(
        description="Reproduce the MATLAB plot_bfield_on_axis_MPEX.m plot in Python."
    )
    parser.add_argument(
        "--config",
        nargs="+",
        default=defaults["config_name"],
        help="Configuration names, for example D3-6 D2-2 D1-1",
    )
    parser.add_argument(
        "--field-models",
        nargs="+",
        default=defaults["field_models"],
        choices=("jackson", "hybrid"),
        help="Field models to compare",
    )
    parser.add_argument(
        "--primary-model",
        default=defaults["primary_model"],
        choices=("jackson", "hybrid"),
        help="Model drawn with a solid line",
    )
    parser.add_argument("--hybrid-nturns", type=int, default=defaults["hybrid_nturns"])
    parser.add_argument("--hybrid-nlayers", type=int, default=defaults["hybrid_nlayers"])
    parser.add_argument("--n-test", type=int, default=defaults["n_test"])
    parser.add_argument("--z-min", type=float, default=defaults["z_limits"][0])
    parser.add_argument("--z-max", type=float, default=defaults["z_limits"][1])
    parser.add_argument("--plot-coils", action="store_true")
    parser.add_argument("--quiet", action="store_true", help="Suppress verbose coil-current output")
    parser.add_argument("--save", type=Path, default=None, help="Save the main figure to this path")
    parser.add_argument(
        "--save-coils",
        type=Path,
        default=None,
        help="Save the optional coil cross-section figure to this path",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Create the figure without calling matplotlib.pyplot.show()",
    )
    return parser


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()

    result = plot_bfield_on_axis_mpex(
        config_name=args.config,
        field_models=args.field_models,
        primary_model=args.primary_model,
        hybrid_nturns=args.hybrid_nturns,
        hybrid_nlayers=args.hybrid_nlayers,
        n_test=args.n_test,
        z_limits=(args.z_min, args.z_max),
        plot_coils=args.plot_coils or args.save_coils is not None,
        verbose=0 if args.quiet else 1,
    )

    if args.save is not None:
        result["fig"].savefig(args.save, dpi=150)
    if args.save_coils is not None and result["coil_fig"] is not None:
        result["coil_fig"].savefig(args.save_coils, dpi=150)

    if not args.no_show:
        import matplotlib.pyplot as plt

        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
