# Validate Find_Orb
An end-to-end testing package that will be able to compare and test orbit prediction accuracy for the Find_Orb open source integrator.
See: https://www.projectpluto.com/find_orb.htm for more Find_Orb related information.

## Installation and Dependencies (for now)
This package currently requires [thor](https://github.com/moeyensj/thor) and *most* of its dependencies and it needs to be running in **development mode**.

Follow these terminal commands to install thor first:

1.
```
git clone https://github.com/moeyensj/thor.git
```
2.
```
cd thor
```
3.
```
git checkout findorb-dev
```
4.
```
python setup.py develop
```
Then check if its properly installed by typing in the python command line,
```
>>> import thor
```
---
To install *Validate Find_Orb* as a python package, cd into the folder once cloned and run this in the terminal,

```
python setup.py develop
```
Then check if its properly installed by typing in the python command line,
```
>>> import validate_findorb
```
---
BIG NOTE:
**This project is still in active development** and will have documentation soon!

If you have any more questions email me here: aidan@b612foundation.org
