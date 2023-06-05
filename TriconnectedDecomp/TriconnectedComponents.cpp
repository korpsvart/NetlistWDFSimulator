/* MEX function developed for integration with the OGDF library
*/

/*
For compilation, run the following command inside this file folder:
mex -DDEBUG -O -largeArrayDims -Iogdf/include/ -Iogdf/build/include/ -Logdf -lOGDF -lCOIN TriconnectedComponents.cpp -v
(remove -DDEBUG flag to avoid debugging messages)

where "ogdf" is the folder containing the library source and the
compiled version (built for example using cmake and then visual studio)
*/


#ifdef DEBUG 
#define DEBUG_MSG(x) x
#else
#define DEBUG_MSG(x)
#endif

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/List.h>
#include <ogdf/decomposition/StaticSPQRTree.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/decomposition/Skeleton.h>

using namespace matlab::data;
using matlab::mex::ArgumentList;


class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {

        ogdf::Graph G;

        buildInputGraph(G, inputs);

        buildTriconnectedComps(G, outputs);

    }

    void buildInputGraph(ogdf::Graph &G, ArgumentList inputs) {
        matlab::data::StructArray inputEdges = inputs[0];
        matlab::data::StringArray inputNodes = inputs[1];

        //Create map for mapping index string name to index number
        //(ogdf library can only identify nodes by their integer index)
        std::map<std::string, ogdf::node> nodesIndexMap;

        int i =0;
        DEBUG_MSG(std::cout << "Printing nodes: \n";)
        for (auto node : inputNodes) {
            std::string nodeName = std::string(inputNodes[i]);
            DEBUG_MSG(std::cout << nodeName << "\n";)
            nodesIndexMap[nodeName] = G.newNode(i);
            i++;
        }

        i=0;
        DEBUG_MSG(std::cout << "Printing edges: \n";)
        for (matlab::data::Struct edge : inputEdges) {
            matlab::data::StringArray endpoints = edge["EndNodes"];
            std::string endpointA = std::string(endpoints[0]);
            std::string endpointB = std::string(endpoints[1]);
            DEBUG_MSG(std::cout << "Edge index " << i << ", (" << endpointA << "," << endpointB << ")\n";)
            G.newEdge(nodesIndexMap[endpointA], nodesIndexMap[endpointB], i);
            i++;
        }


    }

    void buildTriconnectedComps(ogdf::Graph &G, ArgumentList outputs)
    {


        DEBUG_MSG(std::cout <<  "Printing the graph\n";)
        for(ogdf::edge e : G.edges) {
            DEBUG_MSG(std::cout <<"Edge id(??): " << e->index() << ", endpoints: ";)
            DEBUG_MSG(std::cout << e <<"\n";)
        }
        ogdf::Triconnectivity tricComp(G);

        //Create output structure
        const char* field_names[] = {"edges", "type"};
        size_t numberOfComps=tricComp.m_numComp;
        DEBUG_MSG(std::cout << "Number of components: " << numberOfComps << "\n";)

        ArrayFactory factory; //factory for generating output structs and arrays
        StructArray structArray = factory.createStructArray({1, numberOfComps},{"edges", "type"});

        //Structure for returning the endpoints of the edges
        //The number of total edges is M+numberOfComps-1
        //where M is the number of original graph edges
        TypedArray<int16_t> edgesEndpoints = factory.createArray<int16_t>({G.numberOfEdges()+numberOfComps-1, 2});


        //Iterate over the triconnected comps
        int k=0; //index for the components (some components may be empty)
        for(int i=0; i<numberOfComps;i++)
        {
            DEBUG_MSG(std::cout << "component: " << i << "\n";)
            //Get i-th component
            ogdf::Triconnectivity::CompStruct comp = tricComp.m_component[i];
            ogdf::List<ogdf::edge> edgesList = comp.m_edges;
            //Convert edges into some useful representation, and store them
            std::vector<int> edgesIds;
            std::array<ogdf::node,2> endpoints;
            for(ogdf::edge e : edgesList)
            {
                    DEBUG_MSG(std::cout << e->index() << ", endpoints: ";)
                    DEBUG_MSG(std::cout << e << "\n";)
                    edgesIds.push_back(e->index());
                    


                //Store the endpoints inside the edgesEndpoints structArray
                //They are stored by index position for easier access
                endpoints = e->nodes();
                //Maybe find a way to check if it was already assigned before
                //Just for efficiency
                edgesEndpoints[e->index()][0]=endpoints[0]->index();
                edgesEndpoints[e->index()][1]=endpoints[1]->index();

            }
         
            size_t numberOfEdges=edgesIds.size();
            if (numberOfEdges>0) //avoid empty components
            {
                structArray[k]["edges"]=factory.createArray<int>({1, numberOfEdges}, &edgesIds[0],edgesIds.data() + edgesIds.size());
                structArray[k]["type"]=factory.createScalar<int16_t>((int16_t)comp.m_type);
                k++;
            }

        }

        outputs[0] = structArray;
        outputs[1] = factory.createScalar(k); //return also the number of non-empty components
        outputs[2] = edgesEndpoints;

       
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {

    }
};