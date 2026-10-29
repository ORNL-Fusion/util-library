#include "mex.h"
#include <math.h>
#include <stdbool.h>

static const mxArray* get_required_field(const mxArray* s, const char* name) {
    const mxArray* f = mxGetField(s, 0, name);
    if (f == NULL) {
        mexErrMsgIdAndTxt("calc_psi_mex:missingField", "Missing required field '%s'.", name);
    }
    return f;
}

static double get_scalar_field(const mxArray* s, const char* name) {
    const mxArray* f = get_required_field(s, name);
    if (!mxIsDouble(f) || mxIsComplex(f) || mxGetNumberOfElements(f) != 1) {
        mexErrMsgIdAndTxt("calc_psi_mex:badField", "Field '%s' must be a real double scalar.", name);
    }
    return mxGetScalar(f);
}

static const double* get_double_field_ptr(const mxArray* s, const char* name) {
    const mxArray* f = get_required_field(s, name);
    if (!mxIsDouble(f) || mxIsComplex(f)) {
        mexErrMsgIdAndTxt("calc_psi_mex:badField", "Field '%s' must be a real double array.", name);
    }
    return mxGetPr(f);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 3 || nrhs > 4) {
        mexErrMsgIdAndTxt("calc_psi_mex:nrhs", "Usage: [psi,ierr,dpsidr,dpsidz] = calc_psi_mex(g,R,Z,nowarn)");
    }
    if (nlhs > 4) {
        mexErrMsgIdAndTxt("calc_psi_mex:nlhs", "Too many output arguments.");
    }

    const mxArray* g = prhs[0];
    const mxArray* R_in = prhs[1];
    const mxArray* Z_in = prhs[2];

    if (!mxIsStruct(g) || mxGetNumberOfElements(g) != 1) {
        mexErrMsgIdAndTxt("calc_psi_mex:badInput", "First input must be a scalar struct g.");
    }
    if (!mxIsDouble(R_in) || mxIsComplex(R_in) || !mxIsDouble(Z_in) || mxIsComplex(Z_in)) {
        mexErrMsgIdAndTxt("calc_psi_mex:badInput", "R and Z must be real double arrays.");
    }

    mwSize nR = mxGetNumberOfElements(R_in);
    mwSize nZ = mxGetNumberOfElements(Z_in);
    if (nR != nZ) {
        mexErrMsgIdAndTxt("calc_psi_mex:sizeMismatch", "R and Z must have the same number of elements.");
    }

    const double* R = mxGetPr(R_in);
    const double* Z = mxGetPr(Z_in);
    mwSize npts = nR;

    const double* rgrid = get_double_field_ptr(g, "r");
    const double* zgrid = get_double_field_ptr(g, "z");
    double dR = get_scalar_field(g, "dR");
    double dZ = get_scalar_field(g, "dZ");
    double ip_sign = get_scalar_field(g, "ip_sign");
    mwSize mw = (mwSize)get_scalar_field(g, "mw");
    mwSize mh = (mwSize)get_scalar_field(g, "mh");

    const mxArray* coeffs = get_required_field(g, "bicub_coeffs_inv");
    if (!mxIsStruct(coeffs) || mxGetNumberOfElements(coeffs) != 1) {
        mexErrMsgIdAndTxt("calc_psi_mex:badField", "g.bicub_coeffs_inv must be a scalar struct.");
    }

    const double* c00 = get_double_field_ptr(coeffs, "c00");
    const double* c10 = get_double_field_ptr(coeffs, "c10");
    const double* c20 = get_double_field_ptr(coeffs, "c20");
    const double* c30 = get_double_field_ptr(coeffs, "c30");
    const double* c01 = get_double_field_ptr(coeffs, "c01");
    const double* c11 = get_double_field_ptr(coeffs, "c11");
    const double* c21 = get_double_field_ptr(coeffs, "c21");
    const double* c31 = get_double_field_ptr(coeffs, "c31");
    const double* c02 = get_double_field_ptr(coeffs, "c02");
    const double* c12 = get_double_field_ptr(coeffs, "c12");
    const double* c22 = get_double_field_ptr(coeffs, "c22");
    const double* c32 = get_double_field_ptr(coeffs, "c32");
    const double* c03 = get_double_field_ptr(coeffs, "c03");
    const double* c13 = get_double_field_ptr(coeffs, "c13");
    const double* c23 = get_double_field_ptr(coeffs, "c23");
    const double* c33 = get_double_field_ptr(coeffs, "c33");

    plhs[0] = mxCreateDoubleMatrix(1, npts, mxREAL);
    plhs[1] = mxCreateLogicalMatrix(1, npts);

    double* psi = mxGetPr(plhs[0]);
    mxLogical* ierr = mxGetLogicals(plhs[1]);

    double* dpsidr = NULL;
    double* dpsidz = NULL;
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(1, npts, mxREAL);
        plhs[3] = mxCreateDoubleMatrix(1, npts, mxREAL);
        dpsidr = mxGetPr(plhs[2]);
        dpsidz = mxGetPr(plhs[3]);
    }

    const double epsilon = 1e-10;
    double r0 = rgrid[0];
    double z0 = zgrid[0];
    double dR_inv = 1.0 / dR;
    double dZ_inv = 1.0 / dZ;

    for (mwSize i = 0; i < npts; ++i) {
        double Ri = R[i];
        double Zi = Z[i];

        bool bad = mxIsNaN(Ri) || mxIsNaN(Zi);
        double ir_d = floor((Ri - r0 + epsilon) * dR_inv) + 1.0;
        double iz_d = floor((Zi - z0 + epsilon) * dZ_inv) + 1.0;

        if (!bad && (mwSize)ir_d == mw) {
            ir_d = floor((Ri - r0 - epsilon) * dR_inv) + 1.0;
        }
        if (!bad && (mwSize)iz_d == mh) {
            iz_d = floor((Zi - z0 - epsilon) * dZ_inv) + 1.0;
        }

        mwSignedIndex ir = (mwSignedIndex)ir_d;
        mwSignedIndex iz = (mwSignedIndex)iz_d;

        if (bad || ir < 1 || ir > (mwSignedIndex)mw || iz < 1 || iz > (mwSignedIndex)mh) {
            ierr[i] = true;
            psi[i] = mxGetNaN();
            if (dpsidr != NULL) {
                dpsidr[i] = mxGetNaN();
                dpsidz[i] = mxGetNaN();
            }
            continue;
        }

        ierr[i] = false;

        mwSize idx = (mwSize)(ir - 1) + mw * (mwSize)(iz - 1);
        double dr = (Ri - rgrid[ir - 1]) * dR_inv;
        double dz = (Zi - zgrid[iz - 1]) * dZ_inv;

        double drr = dr * dr;
        double drrr = drr * dr;
        double dzz = dz * dz;
        double dzzz = dzz * dz;

        double drz = dr * dz;
        double drzz = dr * dzz;
        double drzzz = dr * dzzz;
        double drrz = drr * dz;
        double drrzz = drr * dzz;
        double drrzzz = drr * dzzz;
        double drrrz = drrr * dz;
        double drrrzz = drrr * dzz;

        psi[i] = ip_sign * (
            c00[idx]       + c10[idx] * dr    + c20[idx] * drr    + c30[idx] * drrr   +
            c01[idx] * dz  + c11[idx] * drz   + c21[idx] * drrz   + c31[idx] * drrrz  +
            c02[idx] * dzz + c12[idx] * drzz  + c22[idx] * drrzz  + c32[idx] * drrrzz +
            c03[idx] * dzzz + c13[idx] * drzzz + c23[idx] * drrzzz + c33[idx] * drrr * dzzz
        );

        if (dpsidr != NULL) {
            dpsidr[i] = ip_sign * (
                c10[idx]        + 2.0 * c20[idx] * dr    + 3.0 * c30[idx] * drr   +
                c11[idx] * dz   + 2.0 * c21[idx] * drz   + 3.0 * c31[idx] * drrz  +
                c12[idx] * dzz  + 2.0 * c22[idx] * drzz  + 3.0 * c32[idx] * drrzz +
                c13[idx] * dzzz + 2.0 * c23[idx] * drzzz + 3.0 * c33[idx] * drrzzz
            ) * dR_inv;

            dpsidz[i] = ip_sign * (
                c01[idx]         + c11[idx] * dr    + c21[idx] * drr    + c31[idx] * drrr   +
                2.0 * c02[idx] * dz  + 2.0 * c12[idx] * drz  + 2.0 * c22[idx] * drrz  + 2.0 * c32[idx] * drrrz  +
                3.0 * c03[idx] * dzz + 3.0 * c13[idx] * drzz + 3.0 * c23[idx] * drrzz + 3.0 * c33[idx] * drrrzz
            ) * dZ_inv;
        }
    }
}
