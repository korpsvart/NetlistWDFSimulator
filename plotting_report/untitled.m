%Plotting the triconnected graph
EdgeTable = table(endNodes, types, ids, values,(0:numel(values)-1)', ...
'VariableNames',{'EndNodes', 'Type', 'Id', 'Value', 'NumericId'});
numEdges = size(E, 1);
%filter virtual
for i=1:numel(T)
    edges = T(i).edges;
    edgesReal = edges(edges<numEdges);
    edgesVirtual = edges(edges>=numEdges);
    localedges = EdgeTable(edgesReal+1, :);
    %Add the virtual edge
    for j=1:numel(edgesVirtual)
        endpointsNodes = N(endpoints(edgesVirtual(j)+1, :)+1);
        vedge(j, :) = {string(endpointsNodes), 'Virtual', 'Virtual', 'Nan', string(edgesVirtual(j))};
    end
 
    localedges = [localedges; vedge];
    G = graph(localedges);
    h = plot(G, "EdgeLabel",  G.Edges.Id+ ' (index:' +G.Edges.NumericId + ')', 'EdgeFontSize', 25, ...
    'NodeFontSize', 30, 'Interpreter', 'latex', 'LineWidth', 2.5, 'EdgeColor','k', 'NodeColor','r', 'EdgeLabelColor', 'b',...
    'NodeLabelColor', 'b');
    if i==2
    h.YData(2) = h.YData(2)+0.5;
    h.YData(4) = h.YData(4)+0.5;
    end
    highlight(h, 'Edges', G.Edges.Type=='Virtual', 'LineStyle','--');
    clear vedge;
end
