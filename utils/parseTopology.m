function [B,Q,G, orderedEdges] = parseTopology(netlistFilename)

parsingResult = strcat(netlistFilename, '_topology.mat');


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

if recompute
   
    % Computing the incidence Matrix A
    A = full(incidence(G));
    dimA = size(A);
    m = dimA(1); %number of nodes
    n = dimA(2); %number of edges (elements)
    plot(G,'EdgeLabel', G.Edges.Id);

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
    save(parsingResult,'B', 'Q','G');
    
end


end


