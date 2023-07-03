function [orderedEdges, Z, S] = all_parsing(netlistFilename, Fs, refEdgeId)

[B,Q,G,orderedEdges] = parseTopology(netlistFilename);

[Z, S] = getZS(netlistFilename,B,Q,orderedEdges,Fs, G,refEdgeId);


end

