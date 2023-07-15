function [G,EdgeTable] = graphFromNetlist(netlistFilename)

%Return graph G and EdgeTable created by parsing
%and LTSpice circuit netlist



%CommentStyle option to ignore SPICE directives
%Range = 2 to ignore first line (specifies Netlist path)

M = readmatrix(netlistFilename, 'OutputType', 'string',...
    'CommentStyle', {'.'}, 'Range', 2);

% Creating circuit graph
endNodes = M(:, 2:3);
types = extractBetween(M(:, 1), 1, 1);
ids = M(:, 1);
values = M(:, 4);
EdgeTable = table(endNodes, types, ids, values,...
'VariableNames',{'EndNodes', 'Type', 'Id', 'Value'});
G = graph(EdgeTable);

end

