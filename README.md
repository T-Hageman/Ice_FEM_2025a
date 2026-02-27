# <img src="Physics/Visualisation/Icon.png" width="150" height="150" align="left"> IceFEM: A finite element code to simulate fracture and viscous/plastic deformation in ice

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Note: This repository serves to preserve the version of code related to the publication: 
> T. Hageman, *Phase-field fracture modelling of ice: From triaxial tests to ice-cliff collapse*, Computer Methods in Applied Mechanics and Engineering, 2026. https://doi.org/10.1016/j.cma.2026.118857


## Overview

IceFEM is a parallel finite element code designed for high-fidelity simulation of ice mechanics, including:
- **Phase-field fracture modeling** for capturing crack propagation and ice failure
- **Viscous and plastic deformation** using advanced constitutive models
- **Thermal effects** and heat transfer in ice systems
- **Large-scale simulations** with MPI parallelization using PETSc

The code is particularly suited for studying ice cliff collapse, crevasse formation, ice shelf dynamics, and triaxial compression tests on ice samples.

![Example Simulation](Documentation/Fig13.png)

## Features

- **Flexible Element Types**: Support for various 2D and 3D element formulations
- **Advanced Material Models**: Viscoplastic Glen's law, elastic-brittle behavior, phase-field fracture
- **Parallel Computing**: MPI-based parallelization for efficient large-scale simulations
- **Multiple Physics**: Coupled thermomechanical, and fracture models
- **HDF5 Output**: Efficient data storage and restart capabilities
- **JSON Input**: Easy-to-configure simulation parameters
- **VTK Visualization**: Visualization while simulations are ongoing, and output files compatible with a range of post-processing options

## Table of Contents

- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Quick Installation (Linux)](#quick-installation-linux)
  - [Manual Installation](#manual-installation)
- [Usage](#usage)
  - [Basic Example](#basic-example)
  - [Running Test Cases](#running-test-cases)
  - [Parallel Execution](#parallel-execution)
- [Project Structure](#project-structure)
- [Configuration](#configuration)
- [Documentation](#documentation)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Installation

### Prerequisites

- **C++ Compiler**
- **MPI**: OpenMPI or MPICH for parallel execution
- **CMake**
- **PETSc**: Portable, Extensible Toolkit for Scientific Computation
- **VTK**: Visualization Toolkit (optional, only required for visualization during runtime)
- **HDF5**: High-performance data storage
- **Eigen3**: Linear algebra library
- **RapidJSON**: JSON parsing library
- **HighFive**: C++ HDF5 wrapper

### Quick Installation (Linux)

An automated installation script is provided for Ubuntu/Debian systems:

```bash
# Clone the repository
git clone https://github.com/T-Hageman/Ice_FEM_2025a.git
cd Ice_FEM_2025a

# Run the installation script (installs all dependencies)
chmod +x Install_ICEFEM.sh
./Install_ICEFEM.sh

# Build the project
cmake .
make -j$(nproc)

# Create a folder for outputs
mkdir Results
```

The installation script will:
1. Install required system packages (compilers, CMake, etc.)
2. Download and compile PETSc with optimized settings
3. Install Eigen3, HDF5, HighFive, and RapidJSON
4. Configure environment variables in `~/.bashrc`

### Manual Installation

If you prefer to install dependencies manually or are using a different operating system:

1. **Install system packages**:
   ```bash
   sudo apt-get update
   sudo apt-get install build-essential cmake gfortran libopenmpi-dev
   ```

2. **Install PETSc**:
   Follow the [PETSc installation guide](https://petsc.org/release/install/)

3. **Set environment variables**:
   ```bash
   export PETSC_DIR=/path/to/petsc
   export PETSC_ARCH=arch-linux-c-opt
   export LIBRARY_DIR=/path/to/libraries
   ```

4. **Install other dependencies** (Eigen3, HDF5, VTK, etc.) in `$LIBRARY_DIR`

5. **Build IceFEM**:
   ```bash
   cmake .
   make -j$(nproc)
   ```

## Usage

### Basic Example

After building, the executable `IceCode` will be in the build directory. Run a simulation using:

```bash
./IceCode path/to/input.json
```

### Running Test Cases

Several test cases are provided in the `TestCases/` directory:

**2D Triaxial Compression Test**:
```bash
./IceCode ./TestCases/TriAxialCompression/TriAxial2D.json
```

![Example Simulation](Documentation/Tri2D.png)

**3D Triaxial Compression Test**:
```bash
./IceCode ./TestCases/TriAxialCompression/TriAxial.json
```

![Example Simulation](Documentation/Tri.png)

**Ice Cliff Collapse Simulation**:
```bash
./IceCode ./TestCases/IceCliffs/IceCliff.json
```

### Parallel Execution

For parallel simulations using MPI:

```bash
mpirun -np 8 ./IceCode ./TestCases/IceCliffs/IceCliff.json
```

Replace `8` with the desired number of processors. 

### Output and Visualization

Results are saved in the `Results/` directory as HDF5 files:
- `results_*.hdf5`: Field data at each output time step
- `mesh_*.hdf5`: Mesh geometry and topology
- `Backup*.hdf5`: Restart files for continuing simulations

Use the included MATLAB scripts, or a custom HDF5 reader to visualize the results.

## Project Structure

```
Ice_FEM_2025a/
├── main.cpp                    # Main entry point
├── CMakeLists.txt              # Build configuration
├── Install_ICEFEM.sh           # Installation script
├── mesh/                       # Mesh handling
│   ├── ElementTypes/           # Element formulations (Quad, Hex, etc.)
│   └── Groups/                 # Element and node groups
├── Physics/                    # Core physics modules
├── Models/                     # Constitutive models and physics
│   ├── LinearElastic/          # Linear elasticity
│   ├── Fracture/               # Phase-field fracture
│   ├── PoroElasticity/         # Porous media mechanics
│   ├── Thermal/                # Heat transfer
│   └── Materials/              # Material definitions
├── Solvers/                    # Linear and nonlinear solvers
├── InputsOutputs/              # I/O handling
├── TestCases/                  # Example simulations
│   ├── TriAxialCompression/    # Triaxial test setup
│   └── IceCliffs/              # Ice cliff collapse
├── utility/                    # Utility functions
├── Documentation/              # LaTeX documentation
└── Results/                    # Simulation output directory
```

## Configuration

Simulations are configured using JSON input files. Key sections include:

- **Mesh**: Geometry definition, element types, boundary groups
- **Physics**: Active models (mechanical, thermal, fracture)
- **Materials**: Material properties (elastic moduli, viscosity, etc.)
- **Boundary Conditions**: Displacements, forces, temperature constraints
- **Time Stepping**: Time integration scheme, step size, duration
- **Output**: Frequency and fields to save



See the `TestCases/` directory for complete examples.

## Documentation

For detailed information about the implementation details, and advanced usage, please refer to:
- **[Full Documentation (PDF)](Documentation/Doc.pdf)**

## Citation

If you use IceFEM in your research or work, please cite the following paper:

```bibtex
@article{Hageman2026,
  title={Phase-field fracture modelling of ice: From triaxial tests to ice-cliff collapse},
  author={Hageman, T.},
  journal={Computer Methods in Applied Mechanics and Engineering},
  year={2026},
  publisher={Elsevier},
  doi={},
  note={In press}
}
```

> T. Hageman, *Phase-field fracture modelling of ice: From triaxial tests to ice-cliff collapse*, Computer Methods in Applied Mechanics and Engineering, 2026. https://doi.org/

---

**Note**: This code is provided for research purposes, and is provided as-is. While care has been taken to verify the simulation results, the author is not responsible for any unintended errors in the code. For detailed theoretical background and validation, please consult the associated publication.