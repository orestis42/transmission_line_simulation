# Transmission Line Visualization

![image](https://github.com/orestis42/transmission_line_visualization/assets/37120208/713d9099-3275-4012-aafa-24313d799add)


## Overview
This script provides a comprehensive analysis and visualization of the magnetic field around two parallel cylindrical conductors (ideal transmission line) using Python. The conductors are aligned along the z-axis with a constant current flowing through them.

## Problem Description

The script addresses the following tasks related to the electromagnetic field analysis:

1. **Vector Potential Calculation**: Compute the magnetic vector potential at any point in the Cartesian coordinate system considering the distances from the point to the conductors.
2. **Graphical Representation**:
   - Plot the vector potential in the xy-plane.
   - Graph normalized vector potential and equipotential lines.
   - Visualize the magnetic induction B field lines using streamlines.

## Installation

To install the necessary dependencies, use [`pipenv`](https://github.com/pypa/pipenv?tab=readme-ov-file#installation). It is recommended to use a virtual environment to manage your dependencies.

Clone the repository.

```bash
git clone git@github.com:orestis42/transmission_line_visualization.git
```

Navigate to the project directory.
```bash
cd transmission_line_visualization
```

Install dependencies using pipenv.

```bash
pipenv install
```

## Usage

```bash
pipenv run python ./src/transmission_line_visualization/transmission_line_visualization.py
```
