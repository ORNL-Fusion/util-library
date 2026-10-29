"""Synthetic checks for Python EFIT flux-surface geometry."""

import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from EFITutils import calc_bfield, compute_flux_surface_geometry, lookup_flux_surface_geometry, postprocess_gfile, postprocess_gfile_geometry
from EFITutils import _geometry_flux_hash


def make_elliptical_gfile(zdim=2.0):
    n = 65
    r_axis, a, b = 2.0, 0.5, 0.7
    r = np.linspace(1, 3, n)
    z = np.linspace(-zdim/2, zdim/2, n)
    rr, zz = np.meshgrid(r, z, indexing='ij')
    g = {
        'mw': n, 'mh': n, 'xdim': 2.0, 'zdim': zdim,
        'rgrid1': 1.0, 'zmid': 0.0, 'rmaxis': r_axis,
        'zmaxis': 0.0, 'ssimag': 0.0, 'ssibry': 1.0,
        'cpasma': [], 'fpol': np.ones(n),
        'psirz': ((rr - r_axis)/a)**2 + (zz/b)**2,
    }
    return postprocess_gfile(g)


class FluxSurfaceGeometryTests(unittest.TestCase):
    def test_fpol_matches_gfile_nodes_and_vacuum_field(self):
        g = make_elliptical_gfile()
        g['fpol'] = 1.0 + 0.3*np.exp(-((g['pn'] - 0.55)/0.06)**2)
        g['bcentr'] = 0.5
        g['rzero'] = 2.0
        postprocess_gfile(g)
        indices = np.array([8, 20, 35, 50])
        radii = 2.0 + 0.5*np.sqrt(g['pn'][indices])
        field = calc_bfield(g, radii, np.zeros_like(radii))
        self.assertFalse(np.any(field['ierr']))
        np.testing.assert_allclose(field['Bt']*radii, g['fpol'][indices], atol=1e-6)
        self.assertNotIn('fpol_coeffs', g)
        self.assertAlmostEqual(calc_bfield(g, 2.55, 0.0)['Bt']*2.55, g['bcentr']*g['rzero'])

    def test_lookup_resolves_curved_x_profile_between_sparse_levels(self):
        g = make_elliptical_gfile()
        postprocess_gfile_geometry(g, refinement=4, psi_n_levels=[0.1, 0.3, 0.6, 0.9])
        psi_n = 0.45
        radius = 2.0 + 0.5*np.sqrt(psi_n)
        lookup = lookup_flux_surface_geometry(g, radius, 0.0)
        self.assertEqual(lookup['iErr'], 0)
        self.assertLess(abs(lookup['x'] - 0.5*(np.sqrt(psi_n) - 1)), 0.001)

    def test_precomputed_profiles_and_lookup(self):
        g = make_elliptical_gfile()
        g = postprocess_gfile_geometry(g, psi_n_levels=[0.1, 0.25, 0.5, 0.75, 0.95])
        table = g['fluxSurfaceGeometry']
        self.assertEqual(table['xDefinition'], 'R_out(psiN,Z=0) - R_out(1,Z=0)')
        self.assertEqual(table['zSampling'], g['zmaxis'])
        lookup = lookup_flux_surface_geometry(g, 2.0 + 0.5*np.sqrt(0.6), 0.0)
        self.assertEqual(lookup['iErr'], 0)
        self.assertLess(abs(lookup['x'] - 0.5*(np.sqrt(0.6) - 1)), 0.01)
        self.assertLess(abs(lookup['volume'] - 2*np.pi**2*2.0*0.5*0.7*0.6), 0.02)
        self.assertAlmostEqual(lookup['transportMetric'], lookup['dVdx']*lookup['gradXSqAvg'])
        outside = lookup_flux_surface_geometry(g, 2.0 + 0.5*np.sqrt(0.99), 0.0)
        self.assertEqual(outside['errorId'], 'OutsideTable')

        shifted = make_elliptical_gfile()
        shifted['zmaxis'] = 0.2
        rr, zz = np.meshgrid(shifted['R'], shifted['Z'], indexing='ij')
        shifted['psirz'] = ((rr - 2.0)/0.5)**2 + ((zz - 0.2)/0.7)**2
        postprocess_gfile(shifted)
        postprocess_gfile_geometry(shifted, refinement=4, psi_n_levels=[0.1, 0.25, 0.5, 0.75, 0.9])
        lookup = lookup_flux_surface_geometry(shifted, 2.0 + 0.5*np.sqrt(0.8), 0.2)
        expected_x = 0.5*(np.sqrt(0.8 - (0.2/0.7)**2) - np.sqrt(1 - (0.2/0.7)**2))
        self.assertEqual(lookup['iErr'], 0)
        self.assertLess(abs(lookup['x'] - expected_x), 0.01)

    def test_geometry_sidecar_reuse(self):
        g = make_elliptical_gfile()
        with TemporaryDirectory() as directory:
            filename = Path(directory)/'gfile'
            filename.write_text('cache provenance')
            source_stat = filename.stat()
            g['filename'] = str(filename)
            g['source_flux_hash'] = _geometry_flux_hash(g)
            g['source_file_stat'] = (source_stat.st_mtime_ns, source_stat.st_size)
            postprocess_gfile_geometry(g, refinement=4, psi_n_levels=[0.1, 0.25, 0.5, 0.75, 0.9])
            self.assertTrue(Path(str(filename) + '.geometry.npz').is_file())
            g.pop('fluxSurfaceGeometry')
            with patch('EFITutils.compute_flux_surface_geometry', side_effect=AssertionError('recomputed')):
                postprocess_gfile_geometry(g, refinement=4, psi_n_levels=[0.1, 0.25, 0.5, 0.75, 0.9])
            self.assertEqual(g['fluxSurfaceGeometry']['nComputed'], 5)

    def test_area_volume_and_refinement(self):
        g = make_elliptical_gfile()
        axis = compute_flux_surface_geometry(g, 2.0, 0.0)
        self.assertEqual(axis['iErr'], 0)
        self.assertEqual(axis['areaRZ'], 0.0)
        self.assertEqual(axis['refinement'], 4)

        for psi_n in (0.25, 0.8, 0.98):
            R = 2.0 + 0.5*np.sqrt(psi_n)
            geo = compute_flux_surface_geometry(g, R, 0.0)
            fine = compute_flux_surface_geometry(g, R, 0.0, refinement=6)
            expected_area = np.pi*0.5*0.7*psi_n
            self.assertEqual(geo['iErr'], 0)
            self.assertEqual(fine['iErr'], 0)
            self.assertLess(abs(geo['areaRZ']/expected_area - 1), 0.01)
            self.assertLess(abs(geo['volume']/(2*np.pi*2.0*expected_area) - 1), 0.01)
            self.assertLess(abs(fine['areaRZ']/geo['areaRZ'] - 1), 0.01)
            self.assertEqual(fine['refinement'], 6)

        tall = make_elliptical_gfile(zdim=14.0)
        geo = compute_flux_surface_geometry(tall, 2.0 + 0.5*np.sqrt(0.8), 0.0)
        self.assertEqual(geo['refinement'], 5)
        self.assertEqual(geo['iErr'], 0)

    def test_arrays_and_independent_status(self):
        g = make_elliptical_gfile()
        r = np.array([[2.0, 2.25], [2.0 + 0.5*np.sqrt(0.8), 2.0 + 0.5*np.sqrt(0.98)]])
        geo = compute_flux_surface_geometry(g, r, 0.0)
        self.assertEqual(geo['psiN'].shape, r.shape)
        self.assertTrue(np.all(geo['iErr'] == 0))
        self.assertTrue(np.all(geo['refinement'] == 4))
        for i in np.ndindex(r.shape):
            scalar = compute_flux_surface_geometry(g, r[i], 0.0, refinement=4)
            self.assertEqual(geo['areaRZ'][i], scalar['areaRZ'])
            np.testing.assert_array_equal(geo['contourR'][i], scalar['contourR'])

        mixed = compute_flux_surface_geometry(g, [2.0, 2.55, 4.0, np.nan, 2.0 + 0.5*np.sqrt(0.8)], 0.0)
        self.assertEqual(mixed['iErr'].tolist(), [0, 1, 1, 1, 0])
        self.assertEqual(mixed['errorId'].tolist(), ['', 'OutsideCore', 'OutsideGrid', 'InvalidPoint', ''])
        self.assertGreater(mixed['psiN'][1], 1)
        self.assertTrue(np.isnan(mixed['areaRZ'][1:4]).all())
        self.assertTrue(all(len(contour) == 0 for contour in mixed['contourR'][1:4]))
        self.assertGreater(mixed['areaRZ'][-1], 0)

    def test_open_and_disconnected_contours(self):
        g = make_elliptical_gfile()
        open_g = dict(g, ssibry=10.0)
        self.assertEqual(compute_flux_surface_geometry(open_g, 2.8, 0.0)['errorId'], 'OpenContour')

        r, z = np.meshgrid(g['R'], g['Z'], indexing='ij')
        island = dict(g)
        island['psirz'] = np.minimum(g['psirz'], ((r - 2.8)/0.15)**2 + (z/0.15)**2)
        postprocess_gfile(island)
        point = 2.8 + 0.15*np.sqrt(0.25)
        self.assertEqual(compute_flux_surface_geometry(island, point, 0.0)['errorId'], 'WrongBranch')
        postprocess_gfile_geometry(island, refinement=4, psi_n_levels=[0.1, 0.25, 0.5, 0.75, 0.9])
        self.assertEqual(lookup_flux_surface_geometry(island, point, 0.0)['errorId'], 'WrongBranch')

    def test_call_wide_validation(self):
        g = make_elliptical_gfile()
        with self.assertRaises(ValueError):
            compute_flux_surface_geometry(g, [2.0, 2.2], np.zeros((2, 1)))
        with self.assertRaises(ValueError):
            compute_flux_surface_geometry(g, 2.2, 0.0, refinement=0)


if __name__ == '__main__':
    unittest.main()
