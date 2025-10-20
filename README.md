# GRASS - Gradient Reflection Aligned Saddle Search

## Usage guide

1. Download the project:
```
git clone https://github.com/TheorChemGroup/GRASS-TS-search.git
```

2. Save your structure into `to_opt.xyz`

3. Write bonds_to_search file (see `./tests` for examples)

4. Run TS search (have to modify tests/da_test, or leave for test calculation):

Use <b>Python 3.12 or 3.10, other versions were not tested</b>

### bonds_to_search file format

The reaction is defined by a configuration file with the following structure:

```
<charge>
<multiplicity>
<solvent>
<DoF1>
<DoF2>
...
<DoFn>
```
* Line 1: charge (int) - The molecular charge
* Line 2: multiplicity (int) - The spin multiplicity  
* Line 3: solvent (string) - The solvent model (e.g., "vacuum", "water")
* Lines 4+: DoF_X - Definitions for Degrees of Freedom (DoFs) that constitute the Reaction Direction Guess (RDG)

**Note:** Any line not starting with b, a, or d is ignored.

#### Defining Degrees of Freedom (DoFs)

Each DoF line specifies a coordinate and its expected direction of change during the reaction.

* Format: [type] [atom_index_1] [atom_index_2] [atom_index_3...] [sign]
* Types:
  * b - Bond (requires 2 atom indices)
  * a - Angle (requires 3 atom indices)
  * d - Dihedral Angle (requires 4 atom indices)
* Sign: The coefficient (+1 or -1) indicates the expected direction of change

#### Example

The following lines define a reaction where a bond between atoms 1 and 2 breaks while a bond between atoms 2 and 3 forms:

```
b 1 2 1
b 3 2 -1
```

This is interpreted as: the (1-2) bond length increases (+1) and the (2-3) bond length decreases (-1).

**Recommendation:** For optimal performance, define DoFs between directly interacting atoms. Using coordinates between distant atoms is discouraged as the associated forces may be weak and poorly correlated with the true reaction coordinate.
# How to use

## From terminal


`TS_find_mirror.py` is usually preferred over `TS_find.py`.

Basic usage:
```
python TS_find.py tests/da_test/to_opt.xyz -tf 0.0001 -tr 100 -p xtb --steps 500 --verbose
```

Parameters:
* -tf - Threshold for force (maximum force along any bond from bonds_to_search must be less than this value)
* -tr - Threshold for relative excess of forces over background (relative excess must be less than this value)
* -p, --program - Quantum chemistry software ("orca" or "xtb")
* --steps, -s - Maximum number of optimization steps
* --verbose - Print detailed output

**Note:** You may use one or both convergence thresholds (-tf and/or -tr).

ORCA-specific example:
```
python TS_find.py tests/da_test/to_opt.xyz -tf 0.0001 -p orca -s 500 -onp 8 -omm 1500 -OPATH /your/path/to/orca -oms "B3LYP ma-def2-TZVP"
```

ORCA-specific parameters:
* -onp, --orca-number-processors - Number of processors for calculation
* -omm, --orca-memory - Memory per processor (MB)
* -OPATH, --ORCA-PATH - Path to ORCA executable (required when -onp > 1)
* -oms, --orca-method-string - Method and basis set string (written after ! in ORCA input file)

## Python usage

from TS_find_mirror import optTS, or from TS_find
```
from TS_find_mirror import optTS #, or from TS_find

optTS(xyz_path: str,
    threshold_force: float=0, 
    threshold_rel: float=0,
    programm=dict(name="xtb"), 
    maxstep:int=7000, 
    print_output:bool=True)
```
where:
- `xyz_path` is path to xyz file 
- `threshold_force=0` is threshold for force (maximum force along any bond from bonds_to_search must be less this value)
- `threshold_rel=0` is threshold for relative excess of forces on  over the background (rel excess must be less this value)
- `maxstep=7000` is maximum number of steps. the search is interrupted if this value is reached
- `print_output` is flag to print output
`threshold_rel` or `threshold_force` must exceed 0. Recommended `threshold_rel`>5, `threshold_force`>0.00002, the less, the more accurately the TS will be found, but at the same time, the longer it will take to find it. If botf `threshold_rel` and `threshold_force` is non-zero both thresholds must converged

