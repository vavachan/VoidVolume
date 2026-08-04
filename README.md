# voidVolume3D

`voidVolume3D` computes the **exact geometric void (free) volume** of a 3-D
packing of possibly-overlapping, possibly-polydisperse spheres. Equivalently,
by subtracting the result from the box volume, it gives the exact volume
occupied by the spheres — with no Monte Carlo sampling or grid discretization
error.

It works by constructing the weighted (radical-plane) Delaunay tessellation
of the sphere centers and integrating, tetrahedron by tetrahedron, the empty
space that is not covered by any sphere.

**Reference algorithm:**
Sastry, S., Corti, D. S., Debenedetti, P. G., & Stillinger, F. H. (1997).
*Statistical geometry of particle packings. I. Algorithm for exact
determination of connectivity, volume, and surface areas of void space in
monodisperse and polydisperse sphere packings.* Physical Review E, 56(5),
5524.

---

## 1. Requirements

- A C++11 compiler (tested with g++; no external libraries required).
- Linux/macOS (the output-directory feature uses POSIX `mkdir`).
- [VMD](https://www.ks.uiuc.edu/Research/vmd/) is optional, needed only if
  you want to *visualize* the Tcl `draw`-command output files described
  below.

## 2. Build

```bash
make
```

or directly:

```bash
g++ voidVolume3D.cpp -std=c++11 -O3 -o voidVolume3D
```

This produces a single executable, `voidVolume3D`.

## 3. Usage

```bash
./voidVolume3D <input_file> <probe_radius> [output_dir] [sample_atom_index]
```

| Argument | Required | Description |
|---|---|---|
| `input_file` | yes | Path to the packing file (format below). |
| `probe_radius` | yes | Non-negative probe radius `r_cut`, added to every particle's radius before the void analysis (see §5). Use `0` for the exact void volume of the packing as given. |
| `output_dir` | no | Directory where all output files are written. Created automatically if it doesn't exist. Defaults to the current directory. |
| `sample_atom_index` | no | 0-based index of one particle whose local Voronoi cell is additionally exported for a quick visual sanity check (see `demo.txt` in §6). Defaults to `0`. Must satisfy `0 <= index < number_of_particles`. |

Example:

```bash
./voidVolume3D N3000GaussianDistribution.dat 0.0 results/
```

## 4. Input file format

A plain-text file with the following layout (this is exactly the format of
the included example, `N3000GaussianDistribution.dat`):

```
<N>
<Lx> <Ly> <Lz>
<x_1> <y_1> <z_1> <r_1>
<x_2> <y_2> <z_2> <r_2>
...
<x_N> <y_N> <z_N> <r_N>
```

- Line 1: `N`, the number of particles.
- Line 2: the three edge lengths of the (orthorhombic) periodic simulation
  box.
- The following `N` lines: each particle's `x y z` center coordinates and
  its radius `r`.

All values may be given in any consistent length unit; the tool does not
rescale them (aside from an internal `SIGMA` hook that is currently fixed
at `1`, kept for possible future unit conversion).

## 5. What `probe_radius` (`r_cut`) means

Before the geometric analysis, every particle's radius is increased by
`probe_radius`:

```
effective_radius = input_radius + probe_radius
```

- `probe_radius = 0` gives the exact void volume of the packing as supplied.
- `probe_radius > 0` dilates every sphere by that amount first. The
  resulting "void" volume is then the space that remains **inaccessible to
  a spherical probe of that radius** rolled through the packing — the
  standard trick used for pore-size / accessible-void analyses (e.g. to
  mimic a solvent or gas molecule of finite size). Sweeping `probe_radius`
  over several runs and collecting the results from
  `r_cut_vs_pore_vol.txt` (append-mode; see §6) builds a
  probe-radius-vs-void-volume curve.

## 6. Output files

All files below are written into `output_dir` (default: current directory).
Files with a VMD Tcl `draw ...` format can be loaded directly in VMD's Tk
console with `source <file>`, or drawn at startup with `vmd -e <file>`.

| File | Format | Contents |
|---|---|---|
| `atoms.txt` | VMD Tcl | `draw sphere` command for every particle, at its dilated radius. |
| `delaunay_edges.txt` | VMD Tcl | Every edge of every Delaunay tetrahedron in the tessellation. **Use this to visually audit the tessellation** — see the PBC caveat in §7. |
| `voidVoronoiVertices.txt` | VMD Tcl | A `draw sphere` marker at every Voronoi vertex (Delaunay tetrahedron circumcenter) classified as void space, after convex-hull trimming. |
| `voronoiEdges.txt` | VMD Tcl | The local Voronoi cell (spheres + edges) of a single particle only (`sample_atom_index`), for a fast, small-scale visual check — rendering this for every particle in VMD would be extremely slow. |
| `demo.txt` | custom CSV | The same single-particle Voronoi-cell geometry as above, in the generic CSV format described below (for use with your own plotting/rendering scripts instead of VMD). |
| `atomsPOVRAY.txt` | custom CSV | Every particle (shape code `1`) **and** the Voronoi (radical-plane) edges of the *entire* system (shape code `2`), in the CSV format below. This is the full-system equivalent of `voronoiEdges.txt`, meant for an external renderer/plotting script rather than VMD. |
| `cav<r_cut>.dat` | VMD Tcl | Void Voronoi vertices, colored and grouped by connected void cluster ("cavity"/pocket), with edges linking neighboring void vertices in the same cluster. One file per `probe_radius` value used (the value is embedded in the filename). |
| `cavities<r_cut>.dat` | custom CSV | The same per-cluster cavity geometry as `cav<r_cut>.dat`, in the CSV format below. |
| `cavityVolume.txt` | TSV | `cluster_index <TAB> void_volume` for every non-empty cavity. |
| `r_cut_vs_pore_vol.txt` | TSV, **append mode** | `probe_radius <TAB> total_void_volume`, one line per run. Because this file is appended to rather than overwritten, running the tool repeatedly with different `probe_radius` values (in the same `output_dir`) accumulates a full probe-radius-vs-void-volume curve. |
| stdout | text | Progress percentage during Delaunay construction, then a final summary: `probe_radius`, box volume, `atoms_volume + void_volume` (should equal the box volume), and the total void volume. |

**Custom CSV format** (`demo.txt`, `atomsPOVRAY.txt`, `cavities<r_cut>.dat`):
each row starts with an integer shape code, followed by that shape's
parameters:

- `1, x, y, z, radius,` — a sphere.
- `2, x1, y1, z1, x2, y2, z2, width,` — a line segment (edge).
- `4, x1, y1, z1, x2, y2, z2, width,` — a thin line segment used for the
  triangulated triangle outlines exported for the sample particle in
  `demo.txt`.

## 7. Periodic boundary conditions — important caveat

Periodic boundary conditions are always applied along all three axes. The
tool does **not** currently check for, or warn about, the following failure
mode inherited from the underlying algorithm:

> If the box is small relative to the particle spacing (or `probe_radius`),
> a particle can end up forming a Delaunay tetrahedron with **another
> particle and that same particle's own periodic image**. This is
> geometrically invalid and will corrupt the void-volume result.

**Always visually check `delaunay_edges.txt` in VMD** after a run,
especially for new/unfamiliar packings or box sizes:

```
vmd -e delaunay_edges.txt
```

Look for any tetrahedron connecting a particle to two entries that both
resolve to the same underlying particle (visually: an edge that appears to
loop a particle back on itself across a box face). If you see this, the box
is too small for the requested `probe_radius`/particle density — use a
larger system (supercell) or reduce `probe_radius`.

A useful diagnostic is the first console line printed after "processing
done": the sum of all Delaunay tetrahedron volumes versus the box volume.
For a geometrically valid tessellation these should match closely; a large
discrepancy is a strong sign of the self-image problem above.

## 8. Performance notes

- The tool performs a breadth-first traversal of the Delaunay neighbor graph
  starting from particle 0, and assumes this reaches every particle (true
  for any ordinary, connected sphere packing). If particles are only
  reachable through in tessellation from disconnected regions, wait — this
  case does not arise in single packed system with periodic boundaries.
- All geometry is computed in `long double` (extended precision) to keep
  the exact-arithmetic predicates (orientation tests, circumcenter
  intersections) numerically robust.
- Runtime scales with the number of particles and their local coordination
  number; the included 3000-particle example runs in well under a minute on
  a modern laptop.

## 9. Known limitations

- Orthorhombic (rectangular) periodic boxes only; no support for triclinic
  cells.
- No automated detection of the particle/self-image tessellation failure
  mode described in §7 — this must currently be checked visually.
- Single-threaded.

## 10. Files in this repository

- `voidVolume3D.cpp` — the tool.
- `Makefile` — `make` to build, `make clean` to remove the binary.
- `N3000GaussianDistribution.dat` — example input: a 3000-particle packing
  with Gaussian-distributed radii.
- `CHANGES.md` — what changed relative to the original script this was
  cleaned up from.
