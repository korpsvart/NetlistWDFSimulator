%Plotting the triconnected graph
EdgeTable = table(endNodes, types, ids, values,(0:numel(values)-1)', ...
'VariableNames',{'EndNodes', 'Type', 'Id', 'Value', 'NumericId'});
numEdges = size(E, 1);
%filter virtual
for i=1:numel(T)
    edges = T(i).edges+1;
    edgesReal = edges(edges<=numEdges);
    edgesVirtual = edges(edges>numEdges);
    localedges = EdgeTable(edgesReal, :);
    %Add the virtual edge
    for j=1:numel(edgesVirtual)
        endpointsNodes = N(endpoints(edgesVirtual(j), :)+1);
        vedge(j, :) = {string(endpointsNodes), 'Virtual', 'Virtual', 'Nan', string(edgesVirtual(j))};
    end
 
    localedges = [localedges; vedge];
    G = graph(localedges);
    h = plot(G, "EdgeLabel",  G.Edges.Id+ ' (index:' +G.Edges.NumericId + ')', 'EdgeFontSize', 20, ...
    'NodeFontSize', 20, 'Interpreter', 'latex', 'LineWidth', 2.5, 'EdgeColor','k', 'NodeColor','r');
    %h.YData(2) = h.YData(2)+0.5;
    %h.YData(4) = h.YData(4)+0.5;
    highlight(h, 'Edges', localedges.Type=='Virtual', 'LineStyle','--');
end
