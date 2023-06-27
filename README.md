
# Netlist Wave Digital Filter Simulator

A MATLAB (R2023a) implementation of a WDF based circuit simulator which can parse a circuit description from an LTSpice netlist.






## Files organization

- Inside the root folder you can find two versions of the main applicative:

    * ```Main.m```
    * ```Main_Tree.m```

    The two versions are described in detail in the following sections.
- The ```utils``` folder contains some MATLAB utility functions, related to the parsing of the topology or the creation of the WDF structure.
- The ```TriconnectedDecomp``` folder contains the code implementing the triconnected decomposition algorithm.
- The ```data``` folder contains some samples netlist files and input/output signals which can be used to test and evaluate the code performance.











## Compilation and running

The MATLAB scripts can simply be opened and executed as they are. However, the triconnected decomposition requires a compiled version of the [OGDF library](https://github.com/ogdf/ogdf/) and of the C++ MEX interface (implemented in ```TriconnectedComponents.cpp```).
An already compiled minimum version of the OGDF library is found in ```TriconnectedDecomp/ogdf/```. Otherwise, if needed, it can be collected from the official repository and compiled (for example using cmake and Visual Studio/XCode or similar).

The MEX interface can instead be compiled by launching this command inside the ```TriconnectedDecomp``` folder:

```
mex -DDEBUG -O -largeArrayDims -Iogdf/include/ -Iogdf/build/include/ -Logdf -lOGDF -lCOIN TriconnectedComponents.cpp -v

```


where ```ogdf``` is the folder containing the library source and the
compiled version.
## Adjustable Netlist Variables

At the start of both scripts there are some variables that can be changed depending on the netlist we want to simulate:


- ```netlistFilename```: specifies the name of the netlist txt file (without the trailing file extension)
- ```refEdgeId```: specifies the (netlist) id of one non-adaptable element in the circuit. If the netlist doesn't have a non-adaptable element, you can specify any element here.
- ```outputPorts```: specifes the list of netlist ids across which you wish to compute the output voltages.
- ```referenceSignalFilenames```: specify the filenames of one or more wave files to be used to validate the results. Leave an empty list if you have no reference file. NB: the wave file must have the same sampling frequency as the input file for the current implementation to work (an additional resampling step could be implemented in the future).



## Main.m (Single junction version)

This version parses an input netlist, creates a computable WDF structure composed of a single adapted junction and simulate the desired output voltages.


## Main_Tree.m (Multiple junctions version)

This version works similarly to the previous one, but the reference circuit is first decomposed into Parallel, Series or Rigid triconnected components. Each of these components becomes a junction in the computable WDF structure.





## "Main" Utility functions (utils)


### ```parseTopology```

- Inputs: 
    - netlistFilename: the netlist filename
- Outputs: 
    - G: the complete circuit graph resulting from the parsing
    - B: fundamental loop matrix
    - Q: fundamental cut-set matrix
    - orderedEdges: a list of edges sorted according to tree-cotree decomposition performed

If the function was already called on a graph which is isomorphic to the current one, it will attempt reusing the previously computed B and Q matrices (which are stored in an appropriate ```.mat``` file).

The result of this function depends only on the topological properties of the network, and not on the circuital elements or values.


### ```getZS```

- Inputs: 
    - netlistFilename: the netlist filename
    - B: fundamental loop matrix
    - Q: fundamental cut-set matrix
    - G: the complete circuit graph resulting from the parsing
    - orderedEdges: a list of edges sorted according to tree-cotree decomposition performed
    - Fs: the sampling frequency
    - refEdgeId: the reference edge id
- Outputs: 
    - Z: the vector of free parameters
    - S: the scattering matrix for a single junction WDF structure


Similarly to ```parseTopology```, it will attempt re-using previously computed Z and S (this is only possible if the network is isomorphic and both the circuital values and the sampling frequency have not been changed).




## "Main_Tree" Utility functions (utils)

### ```getZSFromTriconnected```

- Inputs: 
   - T: contains the triconnected components
   - numEdges: total number of edges in the graph
   - refEdgeIndex: the index of the reference edge
   - E: edges structure containing all the circuits values and informations
   - Fs: sampling frequency
   - endpoints: structure containing the nodes adjacent to each edge
- Outputs: 
    - Tree: an "augmented" version of the input T, containing informations needed for the simulation step 
    - Z: vector of computed free parameters
    - S: a collection of scattering rows for all junctions 

The matrix S contains the scattering rows for each edge of the network. Virtual edges need two scattering rows (one for forward scan, one for backward scan). The forward scan row is stored in the Tree structure field "scatteringUp". Since the number of ports for each junction is not constrained, the matrix S rows should have different number of columns depending on the edge. For simplicity, we allocate a number of columns large enough for the junction with the most ports, and we zero-pad when needed.

This function works mainly as an interface between the main code and the handleComponent recursive function.


### ```handleComponent```

- Inputs: 
    - T: contains the triconnected components
    - N: number of components
    - numEdges: total number of edges in the graph
    - compIndex: index of the current component
    - lastParentEdge: edge used to reach this component
    - E: edges structure containing all the circuits values and informations
    - Z: vector of computed free parameters
    - Fs: sampling frequency
    - depth: depth of the current component
    - endpoints: structure containing the nodes adjacent to each edge
    - overallS: the collection of scattering rows, to be filled during the execution
- Outputs: 
    - T: an "augmented" version of the input T, containing informations needed for the simulation step 
    - Z: vector of free parameters
    - S: a collection of scattering rows for all junctions

Function called by handleComponent, which explores, adapts and compute the scattering matrix S for each component by exploring them in a recursive depth-first fashion.
When first called by getZSFromTriconnected, compIndex should be the index of the component containing the non-adaptable element, and lastParentEdge the edge corresponding to that component.


