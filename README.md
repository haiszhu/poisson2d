# poisson2d

A MATLAB/Fortran solver for the 2D Poisson equation via integral
equations on smooth multiply-connected domains with a composite mesh
for volume potential evaluation.

---

## Structure

```
src/          Fortran 90 source modules
matlab/       MATLAB interface, quadrature utilities, mesh helpers
test/         Driver scripts and unit tests
extern/       Gmsh source tree and prebuilt static library
build/        Compiled objects and libraries (generated)
```

---

## Boundary integral

The BIE discretisation and special quadrature follow the approach of
**Johan Helsing**'s MATLAB demonstration codes for the Laplace equations:

> J. Helsing, demo11b — close evaluation via the RCIP method:
> <https://www.maths.lth.se/na/staff/helsing/Tutor/demo11b.txt>

The Nyström quadrature, panel bookkeeping, and `quadr.m` panel struct
conventions are further aligned with **Alex Barnett**'s BIE2D library,
in particular the close-evaluation panel quadrature:

> A. Barnett, BIE2D — `LapDLP_closepanel.m`:
> <https://github.com/ahbarnett/BIE2D/blob/master/panels/LapDLP_closepanel.m>

Key Fortran modules:

| File | Purpose |
|------|---------|
| `src/utils_mod.f90` | Precision kinds (`r64`, `r128`) |
| `src/geometry_mod.f90` | Panel geometry and normal computation |
| `src/lap2d_mod.f90` | Laplace 2D layer potential kernels |
| `src/specialquad_mod.f90` | Generalized Gaussian / special quadrature |

### GRF disk tests

Two self-contained Fortran drivers verify the BIE solve on the unit disk
against an analytic Green's function representation (GRF):

- `test/test_GRF_disk_r64.f90` — double precision (`r64`)
- `test/test_GRF_disk_r128.f90` — quad precision (`r128`)

Build and run both with:

```bash
make test64 test128
```

Sample output:

```
test_GRF_disk_r64:
  Counts: outside=4652, list1=3712, list2=940
  List1 GRF residual max|u|:   1.4890E-15  log10=  -14.77
  List2 GRF residual max|u|:   3.3594E-15  log10=  -14.45

test_GRF_disk_r128:
  Counts: outside=4652, list1=3712, list2=940
  List1 GRF residual max|u|:   1.8041E-33  log10=  -32.70
  List2 GRF residual max|u|:   3.3153E-33  log10=  -32.45
```

---

## Volume quadrature

Volume potentials require integrating `G(x,y) f(y)` over the domain
interior.  This is handled by a **composite mesh** that partitions the
domain into three region types, each with its own high-order quadrature
scheme.  The code supports (singular/near-singular quadr WIP):

> T. G. Anderson, H. Zhu, S. Veerapaneni,
> *A fast, high-order scheme for evaluating volume potentials on complex
> 2D geometries via area-to-line integral conversion and domain mappings*,
> Journal of Computational Physics, 2023.
> <https://www.sciencedirect.com/science/article/abs/pii/S0021999122007513>

### Three region types

**1. Curved triangles** (near each boundary curve)

Each boundary arc is divided into `nElem` segments.  Each segment spans one
zigzag stride of the shot-point mesh (`compositemesh.m`) and forms a curved
triangle with two straight sides and one curved arc of `Z(t)`.  Volume
quadrature nodes are obtained by pulling Vioreanu–Rokhlin nodes on the
reference triangle through the **Gordon–Hall blended map**, giving
near-machine-precision integration despite the curved edge.

Outputs per boundary (`polygons0`, `polygonsi{k}`):

```
ctrisNodePts   3  × nElem   complex vertices (curve pt, shot pt, next curve pt)
ctrisSidePts   16 × 3 × nElem   GL-sampled points on all three sides
ctrisXq        npols × nElem    x-coordinates of VR quadrature nodes
ctrisYq        npols × nElem    y-coordinates
ctrisWq        npols × nElem    quadrature weights (positive, area element)
```

**2. Straight triangles** (buffer annuli and island buffer zones)

Gmsh (`extern/gmsh`) triangulates the annular buffer region between the
curved-triangle layer and the background grid.  Quadrature is computed by
`trivolquad.m` — affine-mapped Vioreanu–Rokhlin nodes on each Gmsh triangle.

**3. Background boxes** (interior)

A Cartesian grid of axis-aligned boxes covers the bulk of the domain.
`boxvolquad.m` applies tensor-product Chebyshev quadrature (`p × p` nodes
per box) over all active boxes in `BMmask`.

### Assembling quadrature

```matlab
[BMmask, polygons0, polygonsi] = compositemesh(curves, bounds, Nmesh, kress_nodes);
p = 10;

% Curved triangles (all boundaries)
xq_ctri = [polygons0.ctrisXq(:); polygonsi{1}.ctrisXq(:); ...];
wq_ctri = [polygons0.ctrisWq(:); polygonsi{1}.ctrisWq(:); ...];

% Straight triangles (Gmsh + trivolquad)
[xq_tris, yq_tris, wq_tris] = trivolquad(nodes, tris, p);

% Background boxes
[xq_box, yq_box, wq_box] = boxvolquad(BMmask, bounds, gridsize, p);
```

See `test/test_compositemesh.m` for a complete example with three curves
(outer circle + two star-shaped holes) and an area verification check.

### Singular and near-singular quadrature

Singular and near-singular corrections for the volume potential — required
when the target point is on or close to the boundary — will be handled by a
separate modular layer to be integrated into the composite mesh framework.

---

## Build

Requirements: `gfortran-14`, `gcc-14`, `g++-14`, MATLAB (for MEX),
[mwrap](https://github.com/zgimbutas/mwrap).

```bash
make          # core library + MEX + both GRF tests
make core     # Fortran library only (no MEX)
make mex      # MATLAB MEX interface (requires mwrap)
make test64   # double-precision GRF disk test
make test128  # quad-precision GRF disk test
make clean
```

Gmsh static library (one-time build):

```bash
make gmsh
```

---

## License

See `LICENSE`.
