function [tree] = SPQRtree(G, refEdge)


%Step 1:
%Simplify the graph, by merging maximal
%bundles of multiple edges


% We could also do it by using "simplify" function but
% I prefer to have better control
% [G_prime,eind,ecount] = simplify(G);


%The parsing should have already bundled up
%the multiple edges, so we don't need to sort them again

edges = G.Edges;
edges_prime = [];

%Replace bundle with virtual edges, and create components
%from bundled edges

endNodes = edges.EndNodes;

[groups] = findgroups(endNodes(:, 1), endNodes(:, 2));


for i=1:max(groups)
    groupEdges = edges(i==groups, :);
    %Create the virtual edge
    virtualEdge.EndNodes = groupEdges(1, :).EndNodes;
    virtualEdge.Type = "virtual";
    virtualEdge.Id = append("Virtual", num2str(i));
    virtualEdge.Value = "0";
    %Add the virtual edge
    groupEdges = [groupEdges; struct2table(virtualEdge)];
    
    %Create the new component
    components{i, :} = groupEdges;
    
    %Create G' simple graph edge list
    edges_prime = [edges_prime; struct2table(virtualEdge)];
    
    
end

end

