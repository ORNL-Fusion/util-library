# G EQDSK FORMAT

Lao 2/7/97

Source: [G EQDSK FORMAT](https://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf)

## 1. Format

The detailed format is in the Fortran source `.f`. The coordinates are right-handed cylindrical `(R, φ, Z)`. Pressure, the poloidal current function, and the q profile use a uniform flux grid from the magnetic axis to the plasma boundary. Poloidal flux uses a rectangular computation grid; plasma boundary and limiter contours are also supplied.

```fortran
      character*10 case(6)
      dimension psirz(nw,nh),fpol(1),pres(1),ffprim(1),
     . pprime(1),qpsi(1),rbbbs(1),zbbbs(1),
     . rlim(1),zlim(1)
c
      read (neqdsk,2000) (case(i),i=1,6),idum,nw,nh
      read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
      read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
      read (neqdsk,2020) current,simag,xdum,rmaxis,xdum
      read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      read (neqdsk,2020) (fpol(i),i=1,nw)
      read (neqdsk,2020) (pres(i),i=1,nw)
      read (neqdsk,2020) (ffprim(i),i=1,nw)
      read (neqdsk,2020) (pprime(i),i=1,nw)
      read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
      read (neqdsk,2020) (qpsi(i),i=1,nw)
      read (neqdsk,2022) nbbbs,limitr
      read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
c
 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
```

## 2. Variables

| Variable | Meaning and units |
| --- | --- |
| `CASE` | Identification character string. |
| `NW` | Number of horizontal R grid points. |
| `NH` | Number of vertical Z grid points. |
| `RDIM` | Horizontal dimension of the computational box, in meter. |
| `ZDIM` | Vertical dimension of the computational box, in meter. |
| `RLEFT` | Minimum R of the rectangular computational box, in meter. |
| `ZMID` | Z at the center of the computational box, in meter. |
| `RMAXIS` | R of the magnetic axis, in meter. |
| `ZMAXIS` | Z of the magnetic axis, in meter. |
| `SIMAG` | Poloidal flux at the magnetic axis, in Weber /rad. |
| `SIBRY` | Poloidal flux at the plasma boundary, in Weber /rad. |
| `RCENTR` | R at which vacuum toroidal magnetic field `BCENTR` is specified, in meter. |
| `BCENTR` | Vacuum toroidal magnetic field at `RCENTR`, in Tesla. |
| `CURRENT` | Plasma current, in Ampere. |
| `FPOL` | Poloidal current function, `F = RBT`, in m-T, on the flux grid. |
| `PRES` | Plasma pressure, in nt / m2, on the uniform flux grid. |
| `FFPRIM` | `FF’(ψ)`, in (mT)2 / (Weber /rad), on the uniform flux grid. |
| `PPRIME` | `P’(ψ)`, in (nt /m2) / (Weber /rad), on the uniform flux grid. |
| `PSIZR` | Poloidal flux at the rectangular grid points, in Weber / rad. |
| `QPSI` | q values on the uniform flux grid from axis to boundary. |
| `NBBBS` | Number of boundary points. |
| `LIMITR` | Number of limiter points. |
| `RBBBS` | R coordinates of boundary points, in meter. |
| `ZBBBS` | Z coordinates of boundary points, in meter. |
| `RLIM` | R coordinates of the surrounding limiter contour, in meter. |
| `ZLIM` | Z coordinates of the surrounding limiter contour, in meter. |

The source lists `PSIZR` in its variable descriptions but uses `psirz` in the Fortran declaration and read statement above.

## 3. Toroidal Current Density

The source relates toroidal current density to `P’(ψ)` and `FF’(ψ)` as follows:

```text
JT (Amp/m2) = R P’(ψ) + FF’(ψ) / R
```
