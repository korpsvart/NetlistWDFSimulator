function [Tree, Z] = buildTree(T, numEdges, refEdgeIndex, E, Fs)
%Build the SQPR Tree from triconnected components, starting from
%the reference edge

%T contains the triconnected components structs
%numEdges is the total number of edges in the graph
%refEdgeIndex is the index of the reference edge, from which
%we determine the root component of the tree
%E is the edges struct containing all the circuits values and informations
%needed for the adaptation process
%Fs is sampling frequency (for adaptation)

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


 %Port impedances vector (diagonal matrix)
M = max([T.edges]); %find the total number of edges, including virtual ones
Z = zeros(M, 1);

[Tree,Z] = exploreComponent(T, N, numEdges, rootCompIndex, refEdgeIndex, E, Z, Fs, 1);







end