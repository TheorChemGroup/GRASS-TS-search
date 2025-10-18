## General how-to
1. Download the project:
```
git clone https://github.com/ER-Elivian/ts_search.git
```
2. Save your structure into to_opt.xyz

3. Write bonds_to_search file (see ./tests for examples)

4. Run TS search (have to modify tests/da_test, or leave for test calculation):

I use <b>Python 3.12, works om 3.10, earlier versions not tested</b>

### bonds_to_search structure:
```
charge
multiplicity
solvent
DoF1
DoF2
...
DoFn
```
Where charge and multiplicity are integers, solvent may be "vacuum" for no-solvent
Any strings after 3-rd begins not with "b","a" or "d" are inored.

DoF string describes DoF change during reaction.
Use b for bond, a for angle and d for dihedral angle. 
For example, "b 1 2 1" "b 3 2 -1"  means "bond between 1-st and 2-nd athoms and bond 
between 2-nd and 3-rd during the reactions altering with same speed but in opposite direction". If 1-2 breaking, 2-3 forming, if 1-2 forming, 2-3 breaking.

It's not recommended to use bonds between far atoms because such atoms not acting directly and force between that atoms suffring from forces from atoms that are closer to that atoms.
# How to use:
## from terminal
### `TS_find_mirror.py`

```
python TS_find.py tests/da_test/to_opt.xyz -tf 0.0001 -tr 100 -p xtb --steps 500 --verbose
```


where: 
* `-tf` is threshold for force (maximum force along any bond from bonds_to_search must be less this value)
* `-tr` is threshold for relative excess of forces on  over the background (rel excess must be less this value)
* `-p`, `--programm` is the used software ("orca" or "xtb")
* `--steps`, `-s` is maximum number of steps. the search is interrupted if this value is reached
* `--verbose` is flag to print output


```
python TS_find.py tests/da_test/to_opt.xyz -tf 0.0001 -p orca -s 500 -onp 8 -omm 1500 -OPATH /your/path/to/orca -oms "B3LYP ma-def2-TZVP"
```
<b>NOTE:</b> you may use one or both thresholds

where: 
* `-onp`, `--orca-number-processors` is number of procesors using to calculation
* `-omm`, `--orca-memory` is amount of memory per processor  
* `-OPATH`, `--ORCA-PATH` is ORCA `PATH`. It is necessary if `-onp` value is >1
* `-oms`, `--orca-method-string` is string containing method and basis, it will be written in ORCA's input file after `!` sign

### `TS_find.py`
Usually worse than TS_find_mirror, but produces more smooth way and more predictable in changes


## in-code Python using:

from TS_find_mirror import optTS, or from TS_find
```
optTS(xyz_path:str,
    threshold_force:float=0, 
    threshold_rel:float=0,
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

. . .

# `TS_find.ipynb`
This project as .ipynb file. May be non-synchronized with .py version. Here are optimization code and optimization with artificial PES example (example requires only block with `mirror_fn` to execute)
