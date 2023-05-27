function [T, Z] = exploreComponent(T, N, numEdges, compIndex, lastVirtualEdge, E, Z, Fs, depth)

%lastVirtualEdge is the virtual edge from which we came to this
%component (for the root it's the reference edge, so in theory the
%non-linear element)

T(compIndex).depth = depth;

%%Exploration phase
for i=1:numel(T(compIndex).edges)
    
    if T(compIndex).edges(i)==lastVirtualEdge
        continue;
    end

    if T(compIndex).edges(i) >= numEdges %it means the edge is virtual
        %Find the other component with the same virtual edge
        found=false;
        for j=1:N
            if isempty(T(j).children) && isempty(T(j).parent) %to skip already explored nodes
                if any(T(j).edges == T(compIndex).edges(i))
                    %They share the same virtual edge =>
                    %link the components in the tree
                    found = true;
                    T(compIndex).children = [j, T(compIndex).children];
                    T(j).parent = compIndex;
                    %Recursive call
                    [T, Z] = exploreComponent(T, N, numEdges, j, T(compIndex).edges(i), E, Z, Fs, depth+1);
                    break;
                end
            end
        end
        if ~found
            T(compIndex).edges(i)=[];
        end
    end
    

end

%%Adaptation phase
%We explored all the children of this node (or none if it was a leaf)
%so now we can adapt it


%Adapt the real edges
for i=1:numel(T(compIndex).edges)
    if T(compIndex).edges(i) < numEdges %edge is real
        %Get corresponding MATLAB index for the edge
        element_mIndex = T(compIndex).edges(i)+1;
        element = E(element_mIndex);
        Z(element_mIndex) = getZ_element(element, Fs);
    end

end

%Adapt the virtual edge "going up" (lastVirtualEdge)
%or the edge corresponding to the non-linear element
%(it's still lastVirtualEdge if everything was called correctly at the
%root)
element_mIndex = lastVirtualEdge+1;
%Depending on the junction type (0=PARALLEL, 1=SERIES, 2=RIGID)


type = T(compIndex).type;

if type == 0 %PARALLEL
    sum=0;
    for i=1:numel(T(compIndex).edges)
        if T(compIndex).edges(i) ~= lastVirtualEdge %get all the others
            element_mIndex2 = T(compIndex).edges(i)+1;
            sum=sum+(1/Z(element_mIndex2));
        end
    end
    Z(element_mIndex)= 1/sum;
elseif type == 1
    sum=0;
    for i=1:numel(T(compIndex).edges)
        if T(compIndex).edges(i) ~= lastVirtualEdge %get all the others
            element_mIndex2 = T(compIndex).edges(i)+1;
            sum=sum+Z(element_mIndex2);
        end
    end
    Z(element_mIndex)= sum;
elseif type ==2
    %TODO
end






end