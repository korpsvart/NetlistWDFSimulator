function [Tree, Z, S] = getZSFromTriconnected(T, numEdges, refEdgeIndex, E, Fs, endpoints)
%Adapt all junctions and compute scattering
%matrices from triconnected components, starting from
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

refEdgeIndex = refEdgeIndex-1; %decrease index to start from zero
for i=1:N
    if any(T(i).edges==refEdgeIndex)
        rootCompIndex = i;
        break;
    end
end

 %Port impedances vector (diagonal matrix)
M = max([T.edges])+1; %find the total number of edges, including virtual ones
Z = zeros(M, 1);

%"Giant" scattering matrix
%Find highest number of edges for a single comp

%S will contain the rows of the scattering matrix assembled into a single
%big matrix, where the scattering row for element (edge) i
%is found at S(i, :).
%However, virtual edges are shared by 2 components and share the same row.
%So we only store in S the row for the backward scan (i.e. the one
%belonging to the element higher in the tree), while the one for the
%foward scan is stored in the Tree.scatterinUp field.
L = max(arrayfun(@(x) numel(x.edges),T));
S = zeros(M, L);

[Tree,Z, S] = handleComponent(T, numEdges, rootCompIndex, refEdgeIndex, E, Z, Fs, 1, endpoints, S);







end