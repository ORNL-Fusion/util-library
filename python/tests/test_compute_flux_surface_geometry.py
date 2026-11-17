"""Synthetic checks for Python EFIT flux-surface geometry."""

import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from EFITutils import calc_bfield, compute_flux_surface_geometry, lookup_flux_surface_geometry, postprocess_gfile, postprocess_gfile_geometry
from EFITutils import _geometry_flux_hash, calc_flux_surface_labels, lookup_flux_surface_labels


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
        'cpasma': [], 'fpol': np.ones(n), 'rzero': r_axis, 'bcentr': 1/r_axis,
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

    def test_flux_radii_and_metric_transforms(self):
        g = make_elliptical_gfile()
        for p in (0.25, 0.8, 0.98):
            geo = compute_flux_surface_geometry(g, 2+0.5*np.sqrt(p), 0)
            phi = 2*np.pi*(0.7/0.5)*(2-np.sqrt(4-0.25*p))
            radius = np.sqrt(2*phi/np.pi)
            dphi = np.pi*0.5*0.7/np.sqrt(4-0.25*p)
            drdx = 2*dphi*geo['dpsiNdx']/(2*np.pi*radius)
            self.assertEqual(geo['iErr'], 0)
            for name, expected in (('phiTor', phi), ('rTor', radius), ('rSaarelma', radius),
                                   ('dPhiTorDpsiN', dphi), ('drTorDx', drdx)):
                self.assertLess(abs(geo[name]/expected-1), 0.005)
            for suffix, slope, volume in (('RTor', 'drTorDx', 'dVdrTor'),
                                          ('RSaarelma', 'drSaarelmaDx', 'dVdrSaarelma')):
                self.assertAlmostEqual(geo[volume]*geo[slope]/geo['dVdx'], 1)
                self.assertAlmostEqual(geo['grad'+suffix+'SqAvg']/(geo['gradXSqAvg']*geo[slope]**2), 1)
                self.assertAlmostEqual(geo['transportMetric'+suffix]/(geo['transportMetric']*geo[slope]), 1)
        axis = compute_flux_surface_geometry(g, 2, 0)
        self.assertEqual(axis['rTor'], 0)
        self.assertTrue(np.isnan(axis['drTorDx']))

        # Reversing the toroidal field changes sign metadata, not flux magnitudes.
        reversed_g = dict(g, fpol=-g['fpol'], bcentr=-g['bcentr'])
        reverse = compute_flux_surface_geometry(reversed_g, 2+0.5*np.sqrt(0.98), 0)
        self.assertEqual(reverse['toroidalFieldSign'], -1)
        self.assertAlmostEqual(reverse['phiTor'], geo['phiTor'])
        self.assertAlmostEqual(reverse['rTor'], geo['rTor'])

    def test_varying_fpol_and_concave_area_integrals(self):
        from scipy.integrate import quad
        from scipy.interpolate import PchipInterpolator
        from EFITutils import _enclosed_radial_integrals

        # Elliptical area elements integrate analytically in poloidal angle.
        # F=1+0.2*p distinguishes toroidal flux from the geometric Saarelma radius.
        g = make_elliptical_gfile()
        g['fpol'] = 1+0.2*g['pn']
        p = 0.8
        geo = compute_flux_surface_geometry(g, 2+0.5*np.sqrt(p), 0)
        expected = np.pi*0.5*0.7*quad(lambda q: (1+0.2*q)/np.sqrt(4-0.25*q), 0, p)[0]
        self.assertEqual(geo['iErr'], 0)
        self.assertLess(abs(geo['phiTor']/expected-1), 0.005)
        self.assertGreater(geo['rTor'], geo['rSaarelma'])

        # A U-shaped polygon is not star-shaped about the selected axis.
        # Three rectangular pieces give an independent integral of 1/R.
        g['rmaxis'] = 1.25
        contour = np.array([[1, -1], [3, -1], [3, 1], [2.5, 1], [2.5, -0.5],
                            [1.5, -0.5], [1.5, 1], [1, 1], [1, -1]], dtype=float)
        fpol = PchipInterpolator(g['pn'], np.ones_like(g['pn']))
        expected = 0.5*np.log(3)+1.5*np.log(1.5)+1.5*np.log(3/2.5)
        for points in (contour, contour[::-1]):
            phi, integral, error_id = _enclosed_radial_integrals(g, points, fpol)
            self.assertEqual(error_id, '')
            self.assertLess(abs(integral/expected-1), 1e-6)
            self.assertAlmostEqual(phi, integral)

    def test_label_lookup_derivatives_shapes_and_status(self):
        g = make_elliptical_gfile()
        postprocess_gfile_geometry(g, refinement=4, psi_n_levels=[0.1, 0.25, 0.5, 0.75, 0.95])
        p = np.array([[0.25, 0.4], [0.6, 0.75]])
        labels = lookup_flux_surface_labels(g, p)
        spatial = calc_flux_surface_labels(g, 2+0.5*np.sqrt(p), 0)
        self.assertEqual(labels['rTor'].shape, p.shape)
        np.testing.assert_array_equal(labels['iErr'], 0)
        np.testing.assert_array_equal(spatial['iErr'], 0)
        np.testing.assert_allclose(labels['rTor'], spatial['rTor'], atol=1e-8)
        np.testing.assert_allclose(labels['rTor'], labels['rSaarelma'], rtol=1e-12)
        center = lookup_flux_surface_labels(g, 0.6)
        plus = lookup_flux_surface_labels(g, 0.600001)
        minus = lookup_flux_surface_labels(g, 0.599999)
        for radius, slope in (('rTor', 'drTorDx'), ('rSaarelma', 'drSaarelmaDx')):
            finite_difference = (plus[radius]-minus[radius])/(plus['x']-minus['x'])
            self.assertLess(abs(finite_difference/center[slope]-1), 1e-7)
        self.assertLess(abs((plus['phiTor']-minus['phiTor'])/2e-6/center['dPhiTorDpsiN']-1), 1e-7)
        p = np.r_[np.nan, 0, 1, g['fluxSurfaceGeometry']['psiN'][[0, -1]]]
        rejected = lookup_flux_surface_labels(g, p)
        np.testing.assert_array_equal(rejected['iErr'], [1, 0, 0, 1, 1])
        self.assertTrue(np.all(np.isnan(rejected['rTor'][[0, 3, 4]])))
        self.assertEqual(rejected['errorId'].tolist(), ['InvalidPoint', '', '', 'OutsideTable', 'OutsideTable'])
        for key in ('bcentr', 'rzero', 'fpol'):
            changed = dict(g)
            changed[key] = 1.01*g[key]
            with self.assertRaisesRegex(ValueError, 'stale'):
                lookup_flux_surface_labels(changed, 0.6)
        for value in ([1+2j], [], ['0.6']):
            with self.assertRaises(ValueError):
                lookup_flux_surface_labels(g, value)

    def test_separatrix_normalization_and_endpoints(self):
        # Constant F gives equal radii. The enclosed inverse-R integral of
        # an ellipse has an analytic value, including at its outer surface.
        g = make_elliptical_gfile()
        postprocess_gfile_geometry(g, refinement=4, psi_n_levels=[.1, .3, .6, .9])
        table = g['fluxSurfaceGeometry']
        self.assertEqual(table['separatrixIErr'], 0)
        expected_phi = 2*np.pi*.7/.5*(2-np.sqrt(2**2-.5**2))
        self.assertLess(abs(table['phiTorSep']/expected_phi-1), 1e-4)
        self.assertLessEqual(table['separatrixRelativeChange'], 1e-4)
        self.assertGreater(table['rTorSep'], table['rTor'][-1])
        labels = lookup_flux_surface_labels(g, np.array([[0., .6, 1.]]))
        np.testing.assert_array_equal(labels['iErr'], 0)
        np.testing.assert_allclose(labels['rhoTor'][0, [0, 2]], [0, 1])
        np.testing.assert_allclose(labels['rhoTor'], labels['rhoSaarelma'], rtol=1e-12)
        np.testing.assert_allclose(labels['xSaarelma'], labels['rSaarelma']-table['rSaarelmaSep'])
        self.assertEqual(labels['xSaarelma'][0, 2], 0)
        self.assertTrue(np.isnan(labels['drSaarelmaDx'][0, 2]))
        self.assertAlmostEqual(lookup_flux_surface_labels(g, 0)['x'], -.5)
        # Metrics remain interior-only even though the flux-label endpoint exists.
        self.assertEqual(lookup_flux_surface_geometry(g, 2.5, 0)['iErr'], 1)
        # An unavailable edge must not silently normalize to the final sample.
        table.update(phiTorSep=np.nan, rTorSep=np.nan, rSaarelmaSep=np.nan,
                     separatrixIErr=1, separatrixErrorId='SeparatrixNotConverged')
        interior = lookup_flux_surface_labels(g, .6)
        self.assertEqual(interior['iErr'], 0)
        self.assertTrue(np.isnan(interior['rhoTor']))
        self.assertEqual(lookup_flux_surface_labels(g, 1)['errorId'], 'SeparatrixNotConverged')

    def test_close_contours_and_version_upgrade(self):
        from EFITutils import _point_inside_flux_contour

        g = make_elliptical_gfile()
        levels = [0.79, 0.8, 0.800001, 0.81]
        postprocess_gfile_geometry(g, refinement=4, psi_n_levels=levels)
        p = 0.8000005
        radius = 2+0.5*np.sqrt(p)
        table = g['fluxSurfaceGeometry']
        i = np.searchsorted(table['psiN'], p)
        contour = np.column_stack((table['contourR'][i], table['contourZ'][i]))
        self.assertFalse(_point_inside_flux_contour(radius, 0, contour))
        geo = lookup_flux_surface_geometry(g, radius, 0)
        labels = lookup_flux_surface_labels(g, p)
        self.assertEqual(geo['iErr'], 0)
        self.assertAlmostEqual(geo['rTor'], labels['rTor'])
        table['version'] = 1
        with self.assertRaisesRegex(ValueError, 'rebuild'):
            lookup_flux_surface_labels(g, p)
        postprocess_gfile_geometry(g, refinement=4, psi_n_levels=levels)
        self.assertEqual(g['fluxSurfaceGeometry']['version'], 3)

    def test_call_wide_validation(self):
        g = make_elliptical_gfile()
        with self.assertRaises(ValueError):
            compute_flux_surface_geometry(g, [2.0, 2.2], np.zeros((2, 1)))
        with self.assertRaises(ValueError):
            compute_flux_surface_geometry(g, 2.2, 0.0, refinement=0)


if __name__ == '__main__':
    unittest.main()
