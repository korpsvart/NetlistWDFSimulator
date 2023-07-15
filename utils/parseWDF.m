function [orderedEdges, Z, S] = parseWDF(netlistFilename, Fs, refEdgeId, makeGeneratorsReal)

[B,Q,G,orderedEdges] = getTopology(netlistFilename);

if makeGeneratorsReal
    orderedEdges(orderedEdges(:, 1)=="V", 1) = "Vreal";
end

[Z, S] = getZS(netlistFilename,B,Q,orderedEdges,Fs, G,refEdgeId);


end

