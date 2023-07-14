function [B,Q,G, orderedEdges] = getTopology(netlistFilename)

parsingResult = strcat(netlistFilename, '_topology.mat');


% Parsing from LTSpice netlist

[G, ~] = graphFromNetlist(strcat('data/netlist/', netlistFilename , '.txt'));


recompute = true;

if (isfile(parsingResult))
    %Load parsing result if available
    previousData = load(parsingResult, 'B', 'Q', 'G');


    %Check if the graphs are isomorphic
    %If they are, possibly reorder the nodes
    order = isomorphism(G, previousData.G);
    if numel(order)>0 %isomorphic
        G = reordernodes(G,order);
        tree = minspantree(G);

        idx = ismember(G.Edges.Id, tree.Edges.Id);
        idx_n = find(idx);
        cotree = rmedge(G, idx_n);
        orderedEdges = [cotree.Edges.Variables; tree.Edges.Variables];
        
        %Then we can re-use previous B and Q
        B = previousData.B;
        Q = previousData.Q;
        
        recompute = false;
    end
end

if recompute %recompute
   
    [B, Q, orderedEdges] = getBQ(G);
    
    save(parsingResult,'B', 'Q','G');
    
end


end


