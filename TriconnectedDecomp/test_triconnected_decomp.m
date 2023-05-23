%%Code for testing the triconnected components decomposition interface


E = table2struct(G.Edges);
N = table2struct(G.Nodes);
for i=1:numel(E)
    E(i).EndNodes = convertCharsToStrings(E(i).EndNodes);
end
N = string(struct2cell(N));
[T, numComps] = SQPRTree(E, N);
T = T(1:numComps);
numEdges = size(E, 1);
refEdgeIndex = 0; %chosen arbitrarily here, in practice you should give the edge of the non-linear element

%Also pass the element values for adaptation

[Tree, Z] = buildTree(T, numEdges, refEdgeIndex, E, Fs);
