import numpy as np
import cube

# Cube file paths
file_cube_ref = "source/h2_esp.cube"
file_cube_model = "source/h2_esp.cube.mdcm.cube"

# Read cube files
cube_ref, meta_ref = cube.read_cube(file_cube_ref)
cube_model, meta_model = cube.read_cube(file_cube_model)

# Extract data
n = cube_ref.shape
o = np.array(list(meta_ref['org'])[:3])
v = np.array(
    (
        np.array(list(meta_ref['xvec'])),
        np.array(list(meta_ref['yvec'])),
        np.array(list(meta_ref['zvec']))
    )
)
nzero = np.array([
    abs(np.round(o[i]/v[i, i]).astype(int)) for i in range(3)])

# Get (x,z) plane through H2
cube_ref_2d = cube_ref[:, nzero[1], : ]

# Get spatial meshgrid
xarr = np.linspace(o[0], 

