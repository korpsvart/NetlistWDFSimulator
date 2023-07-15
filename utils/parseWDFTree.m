function [Tree,E, Z, S, numComps, outputPorts] = parseWDFTree(netlistFilename, refEdgeId, numOutputs, outputPortsIds, Fs , makeGeneratorsReal)


%% Parsing netlist and building graph G


% Parsing from LTSpice netlist


[G, EdgeTable] = graphFromNetlist(strcat('data/netlist/', netlistFilename , '.txt'));

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

if (makeGeneratorsReal)
    idx = [E.Type]=="V";
    for i=find(idx)
        E(i).Type = "Vreal";
    end
end


%% Adaptation phase
numEdges = size(E, 1);

[Tree, Z, S] = getZSFromTriconnected(T, numEdges, refEdgeIndex, E, Fs, endpoints);



end

