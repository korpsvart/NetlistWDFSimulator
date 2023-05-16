/* MyMEXFunction
    For testing try to call in MATLAB
    SQPRTree(); (input not used now)
*/

/*
For compilation:
mex -O -largeArrayDims -Iogdf/include/ -Iogdf/build/include/ -Logdf -lOGDF -lCOIN SQPRTree.cpp -v

where "tric_matlab" is the folder containing the library source and the
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

        //Nodes
        ogdf::List<std::pair<int,int>> edges{std::make_pair(1,2)};

        //ogdf::customGraph(G,n,edges);


        //ogdf::GraphIO::write(G, "output-manual.svg", GraphIO::drawSVG);
        testTriconnectedComps(outputs);



    }

    void testTriconnectedComps(ArgumentList outputs)
    {

        ogdf::Graph G;
        int n = 2;
        ogdf::randomBiconnectedGraph(G, 4, 10);

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
            structArray[i]["edges"]=factory.createArray<int>({1, numberOfEdges}, &edgesIds[0],edgesIds.data() + edgesIds.size());
            structArray[i]["type"]=factory.createScalar<int16_t>((int16_t)comp.m_type);
        }

        outputs[0] = structArray;

       
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