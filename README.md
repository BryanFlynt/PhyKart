# PhyKart
Processing of MyChron Data for Kart Physics

## Requirements
The requirments are pretty standard with the addition of the pyproj library for [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System) coordinate tranformations.
- Numpy
- SciPy
- Matplotlib
- [Pyproj](https://pypi.org/project/pyproj/) - For cartographic projections and coordinate transformations

## To Run
At the moment you only need to run a single file to generate the plots.
```
> cd src
> python main.py
```

## Other Stuff
There are other files with routines to generate G-G-Velocity envelopes, etc. but they are not used at the moment.
