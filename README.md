# Netlist Wave Digital Filter Simulator

The Netlist Wave Digital Filter (WDF) Simulator is a MATLAB (R2023a) implementation of a circuit simulator based on the Wave Digital Filter method. It can parse circuit descriptions from LTSpice netlist files and simulate the behavior of the circuits.

## File Organization

The project files are organized as follows:

- `Main.m` and `MainTree.m`: These are the main MATLAB scripts for the application, implementing the two different versions of the program (detailed later).
- `utils` folder: Contains MATLAB utility functions related to parsing the circuit topology and creating/simulating the WDF structure. It is further subdivided into
  - `common`: contains utility functions common to both code versions
  - `main`: contains utility functions for the `Main` version
  - `mainTree`: contains utility functions for the `MainTree` version
- `TriconnectedDecomp` folder: Contains the code implementing the triconnected decomposition algorithm.
- `data` folder: Contains sample netlist files and input/output signals for testing and evaluation.

## Compilation and Running

The MATLAB scripts can be executed as they are. However, the triconnected decomposition requires a compiled version of the [OGDF library](https://github.com/ogdf/ogdf/) and of the C++ MEX interface (implemented in `TriconnectedComponents.cpp`).

An already compiled minimum version of the OGDF library is found in `TriconnectedDecomp/ogdf/`. Otherwise, if needed, it can be collected from the official repository and compiled (for example using CMake and Visual Studio/XCode or similar tools).

The MEX interface can instead be compiled by launching this command inside the `TriconnectedDecomp` folder:
```
mex -DDEBUG -O -largeArrayDims -Iogdf/include/ -Iogdf/build/include/
  -Logdf -lOGDF -lCOIN TriconnectedComponents.cpp -v
```

where `ogdf` is the folder containing the library source and the compiled version.

## Adjustable Netlist Variables

At the beginning of both `Main.m` and `MainTree.m`, there are some variables that can be adjusted based on the netlist to be simulated. These variables include:

- `makeGeneratorsReal`: If set to `true`, it will add a negligible resistance value to ideal voltage or current generators to make them adaptable.
- `netlistFilename`: Specifies the name of the netlist file (without the file extension).
- `refEdgeId`: Specifies the (netlist) ID of a non-adaptable element in the circuit.
- `outputPortsIds`: Specifies a list of netlist IDs for computing the output voltages.
- `referenceSignalFilenames`: Specifies the filenames of wave files used for result validation. Leave the list empty if no reference file is available.


## Help and Documentation

For detailed information on the individual functions, please refer to the [documentation](./documentation/Documentation.pdf).

## Authors

- Riccardo Di Bella (riccardo1.dibella@mail.polimi.it)
- Giuseppe Risitano (giuseppe.risitano@mail.polimi.it)


