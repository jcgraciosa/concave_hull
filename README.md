# Concave hull generator
## Description
This package can be used to create a concave hull from a set of points. These points can be in 2D or in 3D (not yet tested).
Reference paper: FIXME - place title and link of paper here

## Libraries/ Requirements
- Python 3.6
- numpy

## Installation Guide
1. `pip install -r requirements.txt`
2. `cd clone_dir`
3. `pip install ./concave_hull`

## Main Functions

1. compute for the concave hull (2D case - tested)
```
FIXME: what should be the main functions?

```

## Scripts and Notebooks

1. concave_test
```
Demonstration of the package's capabilities
```

## Sample Usage

1. Computing the concave hull of a set of 2D points
```python
import concave_hull
"""
Generate set of data here. The package accepts a numpy array containing the data.

2D case - [nrow, 2] data format
3D case - [nrow, 3] data format

See concave_test notebook for sample toy data
"""
save_dir = 'dir_you_want_to_save_plots_for_debugging'
test = concave_hull.ConcHuller(1, data_2d, debug = True, savedir = save_dir)
# note that you only set debug to true and savedir to something valid if you want to debug (duhh)
test.compute_hull()

```
2. Computing the concave hull of a set of 3D points
```python
import concave_hull
"""
Generate set of data here. The package accepts a numpy array containing the data.

2D case - [nrow, 2] data format
3D case - [nrow, 3] data format

See concave_test notebook for sample toy data
"""
save_dir = 'dir_you_want_to_save_plots_for_debugging'
test = concave_hull.ConcHuller(1, data_3d)
# debug option is not yet available for 3D
test.compute_hull()
```

## Test data
Please see concave_test notebook
