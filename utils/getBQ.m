function [B,Q, orderedEdges] = getBQ(G)

    % Computing At and Ac 
    tree = minspantree(G);
    At = full(incidence(tree));
    q = size(At, 2); % q is the number of independent voltages

    %Extract cotree

    %Find edges in the original graph which are part of cotree
    idx = ismember(G.Edges.Id, tree.Edges.Id);

    %Convert logical array into sequential indexes
    idx_n = find(idx);

    %remove edges from original
    cotree = rmedge(G, idx_n);
    Ac = full(incidence(cotree));
    p = size(Ac, 2); % p is the number of independent currents

    % Computing F
    F = pinv(At)*Ac; %F has size q x p 

    % Computing B and Q
    B = [eye(p) -F'];
    Q = [F eye(q)];

    orderedEdges = [cotree.Edges.Variables; tree.Edges.Variables];
end

