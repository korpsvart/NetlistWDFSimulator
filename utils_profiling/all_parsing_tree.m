function [Tree,E, Z, S, numComps, outputPorts] = all_parsing_tree(netlistFilename, refEdgeId, numOutputs, outputPortsIds, Fs)




%% Parsing netlist and building graph G


% Parsing from LTSpice netlist

%CommentStyle option to ignore SPICE directives
%Range = 2 to ignore first line (specifies Netlist path)

M = readmatrix(strcat('data/netlist/', netlistFilename , '.txt'), 'OutputType', 'string',...
    'CommentStyle', {'.'}, 'Range', 2);

% Creating circuit graph
endNodes = M(:, 2:3);
types = extractBetween(M(:, 1), 1, 1);
ids = M(:, 1);
values = M(:, 4);
EdgeTable = table(endNodes, types, ids, values,...
'VariableNames',{'EndNodes', 'Type', 'Id', 'Value'});
G = graph(EdgeTable);

refEdgeIndex = find(EdgeTable.Id==refEdgeId);
outputPorts = zeros(1, numOutputs);
for i=1:numOutputs
    outputPorts(i) = find(EdgeTable.Id==outputPortsIds(i));
end


%% Triconnected components decomposition


%Prepare data for triconnected component decomposition
E = table2struct(EdgeTable);
N = table2struct(G.Nodes);
for i=1:numel(E)
    E(i).EndNodes = convertCharsToStrings(E(i).EndNodes);
end
N = string(struct2cell(N));
[T, numComps, endpoints] = TriconnectedComponents(E, N); %Call C++ interface
T = T(1:numComps);


%% Adaptation phase
numEdges = size(E, 1);

[Tree, Z, S] = getZSFromTriconnected(T, numEdges, refEdgeIndex, E, Fs, endpoints);



end

