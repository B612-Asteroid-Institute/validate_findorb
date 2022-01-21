# Validate Find_Orb
An end-to-end testing package that will be able to compare and test orbit prediction accuracy for the Find_Orb open source integrator.
See: https://www.projectpluto.com/find_orb.htm for more Find_Orb related information.

## Installation and Dependencies

To install *Validate Find_Orb* as a python package, clone this repository, and `cd` into the folder. Then follow the directions below for your use case,

### Installing Dependencies with conda
To install its dependencies in a new `conda` environment, use this command,
```
conda create -n myenv --file requirements.txt python=3.9
```
To install on a preexisting conda environment,
first activate your environment, `conda activate myenv`, then,
```
conda install --file requirements.txt
```

Once all the dependencies have been installed, run this command to install the development version of this package,
```
python setup.py develop --no-deps
```
### Installing Dependencies with pip
Just type in the command line,
```
python setup.py develop
```
**Note:** It is not recommended to use `python setup.py install` for installation since this project is still in development and may have frequent updates.

---
Then check if its properly installed by typing in the python command line,
```
>>> import validate_findorb
```
---
BIG NOTE:
**This project is still in active development** and will have documentation soon!

If you have any more questions email me here: aidan@b612foundation.org
