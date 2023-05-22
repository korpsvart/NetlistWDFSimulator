function [T] = exploreComponent(T, N, numEdges, compIndex)



for i=1:numel(T(compIndex).edges)

    if T(compIndex).edges(i) >= numEdges %it means the edge is virtual
        %Find the other component with the same virtual edge
        for j=1:N
            if isempty(T(j).children) && isempty(T(j).parent) %to skip already explored nodes
                if any(T(j).edges == T(compIndex).edges(i))
                    %They share the same virtual edge =>
                    %link the components in the tree
                    T(compIndex).children = [j, T(compIndex).children];
                    T(j).parent = compIndex;
                    %Recursive call
                    T = exploreComponent(T, N, numEdges, j);
                    break;
                end
            end
        end
    end
    

end


end