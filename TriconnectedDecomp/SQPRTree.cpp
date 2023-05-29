/* MEX function developed for integration with the OGDF library
*/

/*
For compilation:
mex -O -largeArrayDims -Iogdf/include/ -Iogdf/build/include/ -Logdf -lOGDF -lCOIN SQPRTree.cpp -v

where "ogdf" is the folder containing the library source and the
compiled version (built for example using cmake and then visual studio)
*/

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
        std::cout <<"test\n";
        //Actual SQPR Tree code

        ogdf::Graph G;


        buildInputGraph(G, inputs);

        //ogdf::GraphIO::write(G, "output-manual.svg", GraphIO::drawSVG);
        buildTriconnectedComps(G, outputs);



    }

    void buildInputGraph(ogdf::Graph &G, ArgumentList inputs) {
        matlab::data::StructArray inputEdges = inputs[0];
        matlab::data::StringArray inputNodes = inputs[1];
        //int n = matlab::data::ArrayDimensions::getNumElement(inputNodes.getDimensions());


        //Create map for mapping index string name to index number
        //(ogdf library can only identify nodes by their integer index)
        std::map<std::string, ogdf::node> nodesIndexMap;

        //ogdf::node *nodes = (ogdf::node *) malloc( * sizeof(ogdf::node));
        int i =0;
        std::cout << "Printing nodes: \n";
        for (auto node : inputNodes) {
            //nodes[i] = G.newNode(i);
            //nodesIndexMap[inputNodes[i]] = i;
            std::string nodeName = std::string(inputNodes[i]);
            std::cout << nodeName << "\n";
            nodesIndexMap[nodeName] = G.newNode(i);
            i++;
        }

        i=0;
        std::cout << "Printing edges: \n";
        for (matlab::data::Struct edge : inputEdges) {
            matlab::data::StringArray endpoints = edge["EndNodes"];
            std::string endpointA = std::string(endpoints[0]);
            std::string endpointB = std::string(endpoints[1]);
            std::cout << "Edge index " << i << ", (" << endpointA << "," << endpointB << ")\n";
            G.newEdge(nodesIndexMap[endpointA], nodesIndexMap[endpointB], i);
            i++;
        }


    }

    void buildTriconnectedComps(ogdf::Graph &G, ArgumentList outputs)
    {


        std::cout <<  "Printing the graph\n";
        for(ogdf::edge e : G.edges) {
            std::cout <<"Edge id(??): " << e->index() << ", endpoints: ";
            std::cout << e <<"\n";
        }
        ogdf::Triconnectivity tricComp(G);

        //Create output structure
        const char* field_names[] = {"edges", "type"};
        size_t numberOfComps=tricComp.m_numComp;
        std::cout << "Number of components: " << numberOfComps << "\n";

        ArrayFactory factory; //factory for generating output structs and arrays
        StructArray structArray = factory.createStructArray({1, numberOfComps},{"edges", "type"});

        //Iterate over the triconnected comps
        int k=0; //index for the components (some components may be empty)
        for(int i=0; i<numberOfComps;i++)
        {
            std::cout << "component: " << i << "\n";
            //Get i-th component
            ogdf::Triconnectivity::CompStruct comp = tricComp.m_component[i];
            ogdf::List<ogdf::edge> edgesList = comp.m_edges;
            //Convert edges into some useful representation, and store them
            std::vector<int> edgesIds;
            for(ogdf::edge e : edgesList)
            {
                    std::cout << e->index() << ", endpoints: ";
                    std::cout << e << "\n";
                  edgesIds.push_back(e->index());
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

       
    }


    void buildSQPRTree(ogdf::Graph &G, ogdf::edge refEdge, ArgumentList outputs)
    {


        std::cout <<  "Printing the graph\n";
        for(ogdf::edge e : G.edges) {
            std::cout <<"Edge id(??): " << e->index() << ", endpoints: ";
            std::cout << e <<"\n";
        }
        
        
        
        ogdf::StaticSPQRTree sqprTree(G);

        //Get root node (node containing component which has the ref edge)
        ogdf::node rootNode = sqprTree.rootNode();

        //Create output structure
        const char* field_names[] = {"edges", "type", "children"};
        size_t numberOfComps = 5; //adjust this


        ArrayFactory factory; //factory for generating output structs and arrays
        StructArray structArray = factory.createStructArray({1, numberOfComps},{"edges", "type", "children"});


        //get outgoing edges from root node
        std::vector<ogdf::edge> outEdges;
        //rootNode->outEdges(outEdges);

        ogdf::internal::GraphObjectContainer<ogdf::AdjElement> adjEntries = rootNode->adjEntries;



    }



    void testSPQRTree()
    {

        ogdf::Graph G;
        int n = 2;
        ogdf::randomBiconnectedGraph(G, 10, 5);

        ogdf::StaticSPQRTree sqprTree(G);

        //Get root node
        ogdf::node rootNode = sqprTree.rootNode();

        //Get the skeleton
        ogdf::Skeleton& rootEdgeSkeleton = sqprTree.skeleton(rootNode);
        //Get corresponding graph
        ogdf::Graph skelGraph = rootEdgeSkeleton.getGraph();
        for(ogdf::edge e : skelGraph.edges) {
            std::cout << e;
        }
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {

    }
};