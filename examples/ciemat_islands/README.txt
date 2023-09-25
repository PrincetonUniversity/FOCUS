# Prepare Boozer coordinates
You can use the `booz2focus` function in `coilpy` (https://github.com/zhucaoxiang/CoilPy) to export a boundary file.

For example, `booz2focus('boozmn_3p132b.00016.nc', ns=68, focus_file='boozer_ns68_m7n6.boundary')

Here, I use `ns=68` because that is the closest flus surface whose iota = m/n = 7/6 (the target value).

# Prepare a target harmonics file
Now, you should prepare a file that containing the target harmonics, for example, "n6m7.harmonics".

# Prepare FOCUS input file
Prepare a FOCUS input file, there are several variables that should be noted.
```
 input_surf =  'boozer_ns68_m7n6.boundary'   ! define the boundary
 case_surface   =        5              ! 0: general VMEC-like format (Rbc, Rbs, Zbc, Zbs); 1: read axis for knots 5: Booze surface
 weight_bharm   =        1.000D+01      ! weight for Bnm harmonic errors
 input_bharm     =   "n6m7.harmonics"
```

# Run FOCUS
You can now run this case to eliminate islands.

# Misc
The main idea is to read Boozer coordinates and the compute the resonant component in Boozer coordinate. It was first used in https://iopscience.iop.org/article/10.1088/1741-4326/ab3a7c.
