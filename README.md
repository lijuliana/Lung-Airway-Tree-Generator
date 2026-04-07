# `new_airways_work.py` — personalized human lung airway tree generation algorithm

This script **generates a full, personalized human lung airway tree** (on the order of **22–23 generations**) from:

1. **Lung boundary geometry** — a 3D voxel mask of the pleural cavity (currently built as a parametric surface in code).
2. **Trachea / initial airway** — starting segment (position, direction, diameter, length) and the first branching plane.

The algorithm recursively divides available lobe volume between daughter branches using space-dividing planes, flow-based diameters (Murray-type scaling), and optional geometric correction rules (`4a`, `6a`, `7a`, `8a`).

Example:

```bash
pip install numpy scipy sympy matplotlib
```

## How to run

From the project directory:

```bash
python new_airways_work.py
```

## What to configure

### Output file

Near the top of the file, set:

```python
savetofile = "test.pkl"
```

After a successful run, this file contains the **final** segment lists (see [Output](#output)).

### Model and rule toggles

At the top of the script you will find tunables such as:

| Area | Examples |
|------|----------|
| Termination / flow | `flow_rate_threshold` |
| Geometry of branches | `initial_length_to_diameter`, `initial_rotation_angle`, `min_distance_to_length`, `max_distance_to_length` |
| Volume splitting | `volume_dividing_threshold`, `increased_volume_dividing_threshold` |
| Grid refinement | `min_grid_dimension`, `max_grid_dimension` |
| Supplementary rules | `enable_4a`, `enable_6a`, `enable_7a`, `enable_8a` |

Diameter exponent `n` is adjusted in `generate_children` (the implementation uses **3.2** in one segment of the calculation and **2.8** elsewhere, matching the in-file comment about current experimental settings).

### Lung boundary (personalized geometry)

The voxel grid `points` is filled in the **nested loops** after `points = auto_dict()` (search for `for x in range(-72, 73):`). By default, voxels satisfy a ** quartic surface** inequality in normalized coordinates; `total_volume` is accumulated alongside.

To use **another** lung shape, replace that loop body with your own test (inside/outside mask). The main project `README.md` includes examples (sphere, cube) for the same `points[x/4][y/4][z/4] = 1` pattern.

### Trachea and first segment (“trunk”)

The root branch is created as `trunk = Branch(...)`. Important fields:

- **Endpoints:** e.g. from `(0, 0, 0)` to `(0, 0, 7.5)` along **+z**
- **`trachea_diameter`:** scales flow-related child diameters (default `1.8`)
- **`grid_scale`:** initial voxel resolution factor for the root (default `4`)
- **`Plane(..., normal_vector=...)`:** defines the initial branching plane

Adjust these to match **your** trachea pose and size in the same coordinate system as `points`.

## Output

When the run **finishes**, `savetofile` is a **pickle** containing:

```python
[segments_x, segments_y, segments_all]
```

- **`segments_x`:** 2D projection for plotting: each entry is `((x_begin, x_end), (z_begin, z_end), diameter)`.
- **`segments_y`:** `((y_begin, y_end), (z_begin, z_end), diameter)`.
- **`segments_all`:** full 3D tube axis: `((x_begin, x_end), (y_begin, y_end), (z_begin, z_end), diameter)` (values rounded for storage).

**Loading in Python:**

```python
import pickle

with open("test.pkl", "rb") as f:
    segments_x, segments_y, segments_all = pickle.load(f)
```

### Visualization

The script calls `plot_branches(segments)` at the end, which opens **matplotlib** figures: two projections with the analytic boundary curve overlaid, and branch segments drawn with linewidth proportional to diameter.

## Periodic saves during a long run

Every 1000 branches the code writes progress to `savetofile`. Because the file is opened twice in succession, the **intermediate** pickle may only contain `[branch_count, fail_flow_rate]`, not the full segment list. Treat the **final** write after `generate_children` completes as the authoritative tree export.

## Related files

See the project **`README.md`** for notes on `airways_graph.py`, `pig_geometries.py`, supplementary rules, and other geometry recipes.
