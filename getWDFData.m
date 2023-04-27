function [matrixBQ] = getWDFData(netlistFilename)

addpath utils


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

% Returning B or Q matrix
if (q <= p)
    fprintf('Returning matrix B\n');
    matrixBQ = [eye(p) -F']; %B
else
    fprintf('Returning matrix Q\n');
    matrixBQ = [F eye(q)]; %Q
end


orderedEdges = [cotree.Edges.Variables; tree.Edges.Variables];

%Print ordered elements info
fprintf('Ordered elements:\n');
disp([cotree.Edges; tree.Edges]);

%Get Z
ZString = getZString(orderedEdges);

fprintf('Z String:\n %s\n', ZString);


end

