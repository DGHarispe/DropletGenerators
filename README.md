# Microfluidic Droplet Generators on Basilisk

## Overview
This repository contains Basilisk (https://basilisk.fr/) code to simulate droplet generation in microfluidic devices with immiscible phases. Two canonical geometries are supported:

- **Coaxial (flow-focusing) devices** (modeled in 2D assuming axial symmetry).
- **T-junction devices** (modeled in full 3D). Rectangular or cylindrical capillaries

The code targets quantitative analysis of interfacial dynamics, droplet breakup, and regime transitions under controlled physical and geometrical parameters.

## Features

- Simulation of immiscible two-phase flows
- Support for **2D axisymmetric** and **3D geometries**
- Geometry definition:
  - Import from **STL files** (on T-junction capillaries)
- Parameter-driven execution via external "configuration" file

## Numerical Method

- **Discretization**: FVM on a structured (adaptive) Cartesian grid (Quadtree/Octree).
- **Flow solver**: Projection method for incompressible Navier–Stokes equations.
- **Pressure–velocity coupling**: Fractional step method with pressure Poisson equation.
- **Interface capturing**: Volume of Fluid (VOF).
- **Curvature / surface tension model**: Continuum Surface Force (CSF) using Height-Functions.  
- **Time integration**: Explicit/semi-explicit scheme with adaptive timestep.
- **Solid definition (if exists)**: Embedded Boundary Method (and signed distance function when using STLs)

**Related Links**:
- https://basilisk.fr/sandbox/Antoonvh/The_Tree-Grid_Structure_in_Basilisk
- https://basilisk.fr/src/navier-stokes/centered.h
- https://basilisk.fr/src/two-phase.h
- https://basilisk.fr/src/vof.h
- https://basilisk.fr/src/tension.h
- https://basilisk.fr/src/embed.h
- https://basilisk.fr/src/distance.h

## Physical Model

The simulations consider:

- Two immiscible fluids
- Interfacial tension effects
- Different viscosities and densities per fluid
- microchannels geometry

## Geometry Handling

### Coaxial Devices (2D)

- Assumes axial symmetry along the main capillary axis
- Reduced computational cost

### T-Junction Devices (3D)

As specified before, two types of capillaries are available, rectangular and cylindrical capillaries.
The solids are "constructed" by using two different approeaches:

#### STL-based geometries
Import pre-defined capillary geometries from STL files

## Installation

### Requirements

- C compiler
- Basilisk: https://basilisk.fr/src/INSTALL
- [Optional] MPI (for parallel execution)

### Build examples

First, follow the instructions in the Basilisk Installation Guide. Then, use the appropriate compile command based on the type of case you want to compile.

### Coaxial Devices (2D)

No MPI:
```
qcc -source -DDISPLAY=1 -DTRACE=2 source.c
gcc -O2 -Wall _source.c -o source -I$BASILISK -L$BASILISK/gl -L$BASILISK/wsServer -lglutils -lfb_tiny -lws -lGLU -lm
```

### T-Junction Devices (3D)

No MPI:
```qcc -source -grid=octree seToTJ.c
gcc -O2 -Wall _seToTJ.c -o seToTJ -I$BASILISK -L$BASILISK/gl -L$BASILISK/wsServer -lglutils -lfb_tiny -lws -lGLU -lm 
```

MPI:
```
qcc -source -D_MPI=1 -grid=octree seToTJ.c
mpicc -O2 -Wall _seToTJ.c -o seToTJ -I$BASILISK -L$BASILISK/gl -L$BASILISK/wsServer -lglutils -lfb_tiny -lws -lGLU -lm 
```

Related:
https://basilisk.fr/Tutorial#minimal-program

## Usage

1. Create a folder named output next to source or seToTJ once compiled. That folder will be used to place the outputs of the simulations.
2. Execute source/seToTJ.

## Parameters (parameters.txt)

The simulations uses the configuration provided via parameters.txt file. This file is read by functions declared in paramsDict.h and takes the values of certain quantities like viscosities, densities, etc.

### Structure

- Plain text file with key–value pairs
- Parsed at runtime and mapped to internal variables
- Comments can be added using the # character at the start of a line.

Expected format:
```
<parameter_name>: <value>
```

Example:
```
mu_inner: 1e-3
```

### Current supported parameters (as defined by paramsDict.h)
```
NCells: Starting domain size in cells (NCells^2).
maxRlevel: Max refinement level.
q1: Flow rate inner phase.
q2: Flow rate outer phase.
SIGMA: Surface tension.
radius_in: radius of inner capillary for coaxial droplet generator.
radius_out: radius of outer capillary for coaxial droplet generator.
cwidth: Capillary channel width for rectangular capillaries of T-Junction.
cheight: Capillary channel height for rectangular capillaries of T-Junction.
uemax: Max discretisation error for the velocity field.
size_domm: Side length of the (square) domain.
rho1_val / rho_in: Density of the inner phase.
rho2_val / rho_out: Density of the outer phase.
mu1_val / mu_in: Viscosity of the inner phase.
mu2_val / mu_out: Viscosity of the outer phase.
oxc: Origin location on X. 
oyc: Origin location on Y.
ozc: Origin location on Z.
DTi: Max time step.
tmax: Max time for the simulation.
csTol: Max discretisation error for the volume fraction field for solids.
fTol: Max discretization error for the volume fraction field representing the phase interface.
```

The max discretisation errors are used to control the refinement in different regions of the domain, acording to different fields.
There are other parameters defined but currently are not being used on the example codes.

## STL usage
1. Create or obtain STL geometry
2. Place next to seToTJ.c or source.c files.

## Simulation Outputs

### HTG Files
HTG (Hyper Tree Grid) files store hierarchical, adaptive Cartesian grid data, commonly used in AMR-based simulations.

The available code will export the following fields:
- Phases associated volume fraction field (f)
- Pressure field (p)
- Components of velocity field (u)
- Solid associated volume fraction field (cs)
- Solid associated area fraction field (fs)

In ParaView, HTG files can be directly loaded for visualization and analysis of the saved fields.

The subroutines used to export the fields to HTG files is based on a slight modification of the codes found here:

http://basilisk.fr/sandbox/sander/output_htg.h
http://basilisk.fr/sandbox/sander/output_pvd.h

### Basilisk snapshots
In Basilisk, simulation snapshots can be saved on "dump files". These include including grid structure, scalar/vector fields, and solver state. They are used to restart simulations or post-process results. Dumps preserve the adaptive mesh and ensure exact reconstruction of the simulation state.

Related:
https://basilisk.fr/src/output.h#dump-basilisk-snapshots

## Limitations

- No contact angle modeling
- Axisymmetric assumption for coaxial devices
- Other?

## Contributing

Contributions are welcome. Open an issue or submit a pull request with a clear description.

## License

https://basilisk.fr/src/COPYING

## Validation
The code was validated and used to perform the studies available at:

Harispe, D. G., & Kler, P. A. (2024). Accurate numerical prototypes of microfluidic droplet generators with open source tools. Computers &amp; Fluids, 281, 106366. https://doi.org/10.1016/j.compfluid.2024.106366

## Citation

If you use this code in academic work, please cite our work:

@article{Harispe2024,
  title = {Accurate numerical prototypes of microfluidic droplet generators with open source tools},
  volume = {281},
  ISSN = {0045-7930},
  url = {http://dx.doi.org/10.1016/j.compfluid.2024.106366},
  DOI = {10.1016/j.compfluid.2024.106366},
  journal = {Computers & Fluids},
  publisher = {Elsevier BV},
  author = {Harispe, David Gabriel and Kler, Pablo A.},
  year = {2024},
  month = aug,
  pages = {106366}
}

