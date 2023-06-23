EdgeTable = table(endNodes, types, ids, values,(0:numel(values)-1)', ...
'VariableNames',{'EndNodes', 'Type', 'Id', 'Value', 'NumericId'});
G = graph(EdgeTable);

props = plot(G, "EdgeLabel", G.Edges.Id+ ' (index:' +G.Edges.NumericId + ')', 'EdgeFontSize', 15, ...
    'NodeFontSize', 15, 'Interpreter', 'latex');
props.XData(5) = props.XData(5)+0.9;