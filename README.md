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

For calculation structule (XMOL .xyz) and bonds_to_search files are needed
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
* Line 2: multiplicity (int) - The spin multiplicity. If "auto", multiplicity = 1+|charge|
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
  * sign (expected change): The coefficient (typically +1 or -1) indicates the expected direction of change

#### Example

The following lines define a reaction in water with zero full charge where a bond between atoms 1 and 2 breaks while a bond between atoms 2 and 3 forms:

```
0
water
auto
b 1 2 1
b 3 2 -1
```

This is interpreted as: the (1-2) bond length increases (+1) and the (2-3) bond length decreases (-1).

# How to run calculation

## From terminal


`TS_find_mirror.py` is usually preferred over `TS_find.py`.

Basic usage:
```
python TS_find_mirror.py tests/da_test/to_opt.xyz -tf 0.0001 -tr 100 -p xtb --steps 500 --verbose
```

Parameters:
* -tf - Threshold for force (maximum force along any bond from bonds_to_search must be less than this value)
* -tr - Threshold for relative excess of forces over background (relative excess must be less than this value)
* -p, --program - Quantum chemistry software ("orca" or "xtb")
* --steps, -s - Maximum number of optimization steps
* --verbose - Print detailed output

**Note:** You may use one or both convergence thresholds (-tf and/or -tr).

You may prefer standard thresholds like rms and max gradient and displacement. These thresholds may be used that way:

```
python TS_find_mirror.py tests/da_test/to_opt.xyz -tm standard -tgmax 0.0002 -tgrms 0.00004 -tdmax 0.00018 -tdrms 0.0001
```

By default method uses pre-optimization with all DoFs fixed (and biased if `-sa` is used). You may turn it off by using `-nopo`, but it's usually bad decision. Methods is working with assumation that all forces in structure casued by it's  

ORCA-specific example:
```
python TS_find_mirror.py tests/da_test/to_opt.xyz -tf 0.0001 -p orca -s 500 -onp 8 -omm 1500 -OPATH /your/path/to/orca -oms "B3LYP ma-def2-TZVP"
```

ORCA-specific parameters:
* -onp, --orca-number-processors - Number of processors for calculation
* -omm, --orca-memory - Memory per processor (MB)
* -OPATH, --ORCA-PATH - Path to ORCA executable (required when -onp > 1)
* -oms, --orca-method-string - Method and basis set string (written after ! in ORCA input file)

Orca lacks of harmonic constraints, therefore, xtb is used to biasing starting geometry. You may turn this off by `-noppo` (with turning off starting geometry biasing)  

All other parametres may be shown by running
```
python TS_find_mirror.py --help
``` 

## Python usage

```
from TS_find_mirror import optTS #, or from TS_find

optTS(xyz_path: str,
    thresholds={"mode":"native",
                "force": 0.001, 
                "relative": 8},
    programm=dict(name="xtb"), 
    maxstep:int=7000, 
    print_output:bool=True)
```
where:
- `xyz_path` is path to xyz file 
- `thresholds["mode"]` is mode of threcholds. May be `"native"` or `"standard"`
- `thresholds["force"]=0.001` is threshold for force (maximum force along any bond from bonds_to_search must be less this value)
- `thresholds["relative"]=8` is threshold for relative excess of forces on  over the background (rel excess must be less this value)
- `maxstep=7000` is maximum number of steps. the search is interrupted if this value is reached
- `print_output` is flag to print output
- `thresholds["relative"]` and `thresholds["force"]` must exceed 0 if used. Recommended `thresholds["relative"]`>5, `thresholds["force"]`>0.00002, the less, the more accurately the TS will be found, but at the same time, the longer it will take to find it. If both `thresholds["relative"]` and `thresholds["force"]` is non-zero both thresholds must converged

## What happens when you run GRASS calculation?

Example:
```
python TS_find_mirror.py tests/da_test/to_opt.xyz -tf 0.0001 -tr 100 -sa 0.1 -p orca -oms "PBE0 def2-SVP" --steps 500 --verbose
```

1) DoFs in bonds_to_search are measured. There are two DoFs, both are bonds, 1-11 is 2.043 and 2-14 is 2.406 Angstroems and in bonds_to_search they has +1 coefficients. 
2) After, this values are altering by step_along (-sa keybord) by +0.1 Angstroem each (coefficient multiplied by -sa value) and molecular structure optimizes with constraints eqal to DoFs values. Strucure altering to new DoF values is casued by optimisation in xtb with harmonic constraints and big force constant. step_along is useful to destabilize structure if starting point is minimum point, usually 0.7 is enough. Grass itself **can't** escape stationary point therefore manual or automatical destabilization in reaction direction is needed. If you have alredy destabilized structure (like in da_test), this optimization in xtb may be skipped with -noppo (no pre-preoptimisation) key. Even when structure won't changes with -sa, this optimisation is used because of low time of xtb calculations and yield of closer to optimized geometry than constructed manually one.
3) After optimization in xtb, runs optimization in orca to exclude influence on gradient of interatomic inceractions other than casued by displacement in reaction direction. GRASS is based on assumation that majority of forces in structure casued only by displacement in reaction direcrion and searches in this direction for statioary point. If (and typically only if) you have optimized with current level of theory structure, you may skip this step with -nopo (no preoptimisation) flag. -nopo also turns off pre-preoptimisation in xtb to avoid massive destabilization of structure by optimizing in different level of theory. Usually, if structure not optimized at the same level of theory as used for TS search, GRASS calculation takes at least 30-50% more time. If used program is xtb, only this pre-optimisation is used (pre-preoptimisation isn't needed).
4) After all these preparations, starts GRASS TS search: on each step gradient is calculated, reflected by mirror_fn, altered by alter_grad to avoid topology violation in big structures and used in one of optimizers selected in start of calculation. Typically no_Adam is the best one, list of optimizers can be found in help.

# Tips and tricks

* Please note that atoms are numbered from 1, as in XTB and Chemcraft for both XTB and ORCA calculations

* Some transition states can be described in different ways. For example, the turning of ammonia into the planar state (1 N, 2 H, 3 H, 4 H) can be described as synchronous increase of all HNH angles:

      a 2 1 3 1 
      a 3 1 4 1 
      a 4 1 2 1 
  
  But it can also be described as an increase of the dihedral angle along the NH edge 
      
      d 2 1 3 4 1 

  And if in the case of ammonia there is no difference, then for structures with different substituents it is generally better to use the smallest number of degrees of freedom that describe the TS well enough.

* Writing a correct bonds_to_search for a conformal transformation reaction is usually more difficult. But if you need to find the TS of the such reaction, it's not recommended to include bonds in the bonds_to_search that are not present in the molecular structure. For example, in the case of ammonia (see above), the bonds_to_search DoFs **shouldn't** look like this:

      b 2 3 1
      b 3 4 1
      b 4 1 1

* For large structures (e.g., proteins), it may be reasonable to fit sign coefficients on a smaller structure, achieving the fastest convergence. Typically, a range of 0.7–1.4 is sufficient.

* If two transition states are located close to each other on the surface, it is possible to achieve specific convergence to one state or the other by varying the coefficients. For example, in the case of an asynchronous Diels - Alder reaction (Figure 2 [here](https://www.ch.ic.ac.uk/rzepa/motm/porphyrins/introDA.html)) there are two bonds between atoms x1 and y1 and between x2 and y2, and the such bonds_to_search DoFs
      
      b x1 y1 1
      b x2 y2 0.7

  typically guides to the reaction in which the x1-y1 bond is formed earlier, while the such bonds_to_search DoFs
      
      b x1 y1 0.7
      b x2 y2 1
  
  forces bond x2-y2 to form before x1-y1.

* For small structures thresholds may be less (usually converges from `-tf 0.00004 -tr 8`), but if they bigger than `-tf 0.001 -tr 100`, resulting energy may be wrong even for large structures.