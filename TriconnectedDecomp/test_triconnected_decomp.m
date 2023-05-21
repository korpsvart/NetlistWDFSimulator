%%Code for testing the triconnected components decomposition interface


E = table2struct(G.Edges)
N = table2struct(G.Nodes)
for i=1:numel(E)
    E(i).EndNodes = convertCharsToStrings(E(i).EndNodes);
end
N = string(struct2cell(N))
T = SQPRTree(E, N)