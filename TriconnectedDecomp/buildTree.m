function [Tree] = buildTree(T, numEdges, refEdgeIndex)
%Build the SQPR Tree from triconnected components, starting from
%the reference edge

%T contains the triconnected components structs
%refEdgeIndex is the index of the reference edge, from which
%numEdges is the total number of edges in the graph
%we determine the root component of the tree

N = numel(T);

%Add children and parent field to all 

for i=1:N
    if any(T(i).edges==refEdgeIndex)
        rootCompIndex = i;
        break;
    end
end


%%Start building child and parent relationships



%Initialize fields (will add them to all structs, as empty lists)
T(rootCompIndex).children = [];
T(rootCompIndex).parent = -1;

Tree = exploreComponent(T, N, numEdges, rootCompIndex);





end