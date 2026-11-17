"""Compare Python coordinates and metrics with assess_flux_surface_coordinates.m.

Pass the same gfile and the MATLAB MAT file containing Results. The diagnostic
reads the equilibrium without writing a geometry cache and checks direct and
cached values at the MATLAB sampling points. Optional --output saves JSON and
a comparison plot. Run from any directory; the adjacent library is selected.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import EFITutils as efit


def assess_flux_surface_coordinates(gfile_name, matlab_results):
    """Return maximum relative discrepancies and the compared profiles."""
    reference = loadmat(matlab_results, simplify_cells=True)['Results']
    g = efit.postprocess_gfile(efit.readg_g3d_simple(str(gfile_name)))
    table_reference = reference['table']
    # Matching the full raw equilibrium is required before assigning small
    # code-to-code differences to quadrature or contour extraction.
    source = table_reference['Source']
    names = {'R': 'r', 'Z': 'z'}
    for name in ('psirz', 'ssimag', 'ssibry', 'rmaxis', 'zmaxis', 'R', 'Z',
                 'ip_sign', 'pn', 'fpol', 'rzero', 'bcentr'):
        np.testing.assert_allclose(g[name], source[names.get(name, name)], rtol=1e-13, atol=1e-13)
    base_refinement = int(table_reference['refinement'])
    efit.postprocess_gfile_geometry(g, refinement=base_refinement,
                                    psi_n_levels=table_reference['requestedPsiNLevels'])
    fields = efit._GEOMETRY_FIELDS
    report = {'gfile': str(gfile_name), 'matlab_results': str(matlab_results),
              'base_refinement': base_refinement, 'comparisons': {}}
    profiles = {}
    omp = efit._outboard_flux_map(g, g['zmaxis'], base_refinement)
    for label, key, factor in (('direct', 'Direct', 1), ('fine', 'Fine', 2), ('cached', 'Cached', None)):
        matlab = reference[key]
        p = np.array([point['psiN'] for point in matlab])
        radii = np.array([efit._outboard_radius_and_gradient(g, value, omp)[0] for value in p])
        if factor is None:
            geo = efit.lookup_flux_surface_geometry(g, radii, g['zmaxis'])
        else:
            geo = efit.compute_flux_surface_geometry(g, radii, g['zmaxis'], factor*base_refinement, z_omp=0)
        if np.any(geo['iErr']):
            raise AssertionError(f'{label}: {geo["errorId"]}')
        comparisons = {}
        for name in fields:
            values = np.array([point[name] for point in matlab])
            error = np.max(abs(geo[name]-values)/np.maximum(abs(values), np.finfo(float).tiny))
            comparisons[name] = float(error)
        report['comparisons'][label] = comparisons
        profiles[label] = (p, geo, matlab)
        print(f'\n{label} Python / MATLAB: maximum relative discrepancies')
        for name, error in comparisons.items():
            print(f'  {name:27s} {error:.6g}')
        # This is agreement between implementations on identical grids,
        # not a bound on contour discretization error near an X-point.
        if max(comparisons.values()) > 1e-3:
            raise AssertionError(f'{label} Python/MATLAB discrepancy exceeds 0.1%')
    return report, profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gfile')
    parser.add_argument('matlab_results')
    parser.add_argument('--output', type=Path, help='Output prefix for JSON and PNG files')
    args = parser.parse_args()
    report, profiles = assess_flux_surface_coordinates(args.gfile, args.matlab_results)
    if args.output:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        args.output.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
        fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
        fields = (('dVdx', 'dVdrTor', 'dVdrSaarelma'),
                  ('gradXSqAvg', 'gradRTorSqAvg', 'gradRSaarelmaSqAvg'),
                  ('transportMetric', 'transportMetricRTor', 'transportMetricRSaarelma'))
        titles = ("V' (m²)", r'$\langle|\nabla\lambda|^2\rangle$', "V' metric (m²)")
        p, geo, matlab = profiles['fine']
        for ax, group, title in zip(axes, fields, titles):
            for name, label, color in zip(group, ('x OMP', 'r tor', 'r Saarelma'), ('C0', 'C1', 'C2')):
                ax.plot(p, geo[name], color=color, label=label)
                ax.plot(p, [point[name] for point in matlab], 'o', color=color, fillstyle='none')
            ax.set(xlabel=r'$\psi_N$', ylabel=title)
            ax.grid(True)
        axes[0].legend()
        fig.suptitle('Python (lines) and MATLAB (markers), identical contour refinement')
        fig.savefig(args.output.with_suffix('.png'), dpi=150)
    print('\nPython/MATLAB coordinate comparison passed')


if __name__ == '__main__':
    main()
