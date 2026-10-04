# Views of Curves and Surfaces

Sage code for computing visual events of curves and surfaces in P^3, and the
singularities (flexes, bitangents, duals) of real plane curves. It accompanies
K. Kohn, B. Sturmfels, M. Trager, *Changing views on curves and surfaces*
([arXiv:1707.01877](https://arxiv.org/abs/1707.01877)).

## Contents

- `Events for surfaces.ipynb`: flecnodal, cusp crossing, tritangent and parabolic
  curves in Gr(1, P^3), the Hurwitz threefold, bitangent and principal
  congruences, and image contours.
- `Events for curves.ipynb`: tangential surface, trisecant and edge curves, Chow
  form and secant congruence of a space curve.
- `Singularities of Plane Curves.ipynb`: real flexes and bitangents, dual curves,
  Plücker and Klein formulas, and a gallery of real quartics.
- `visual_events/`: the code used by the notebooks.
  - `binary_forms.py`: coincident root loci (discriminants), which encode the
    contact of a line with a curve or surface.
  - `grassmannian.py`: Plücker coordinates and invariants of families of lines.
  - `surfaces.py`, `space_curves.py`, `plane_curves.py`: the computations.
  - `plotting.py`: real plane curves, including a "disk" view of the whole
    real projective plane.
- `tests/`: checks against classical formulas and the tables of the paper.

## Setup

Install SageMath with conda (here via [Miniforge](https://github.com/conda-forge/miniforge)):

```
conda create -n sage sage python=3.11 jupyterlab
conda activate sage
jupyter lab
```

To use the notebooks in VS Code, register the Sage kernel once and choose
"SageMath" as the notebook kernel:

```
jupyter kernelspec install --user ~/miniforge3/envs/sage/share/jupyter/kernels/sagemath
```

Run the tests with

```
sage -python tests/test_events.py
sage -python tests/test_plane_curves.py
```
