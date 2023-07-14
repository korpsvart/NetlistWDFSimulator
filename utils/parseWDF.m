function [orderedEdges, Z, S] = parseWDF(netlistFilename, Fs, refEdgeId)

[B,Q,G,orderedEdges] = getTopology(netlistFilename);

[Z, S] = getZS(netlistFilename,B,Q,orderedEdges,Fs, G,refEdgeId);


end

