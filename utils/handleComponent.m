function [T, Z, overallS] = handleComponent(T, numEdges, compIndex, lastParentEdge, E, Z, Fs, depth, endpoints, overallS)

%lastParentEdge is the virtual edge from which we came to this
%component (for the root it's the reference edge, so in theory the
%non-linear element)

T(compIndex).depth = depth;
T(compIndex).parentEdge = lastParentEdge;

%%Exploration phase
numCompEdges=numel(T(compIndex).edges);
for i=1:numCompEdges
    
    if T(compIndex).edges(i)==lastParentEdge
        continue;
    end

    if T(compIndex).edges(i) >= numEdges %it means the edge is virtual
        %Find the other component with the same virtual edge
        found=false;
        for j=1:numel(T)
            if isempty(T(j).parentEdge) %to skip already explored nodes
                if any(T(j).edges == T(compIndex).edges(i))
                    %They share the same virtual edge =>
                    %link the components in the tree
                    found = true;
                    %Recursive call
                    [T, Z, overallS] = handleComponent(T, numEdges, j, T(compIndex).edges(i), E, Z, Fs, depth+1, endpoints, overallS);
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
for i=1:numCompEdges
    %Check if edge is real and if it's not the root non-adaptable element
    if T(compIndex).edges(i) < numEdges && T(compIndex).edges(i)~=lastParentEdge
        %Get corresponding MATLAB index for the edge
        element_mIndex = T(compIndex).edges(i)+1;
        element = E(element_mIndex);
        Z(element_mIndex) = getZElement(element, Fs);
    end

end

%Adapt the virtual edge "going up" (lastParentEdge)
%or the edge corresponding to the non-linear element
%(it's still lastParentEdge if everything was called correctly at the
%root)
element_mIndex = lastParentEdge+1;
%Depending on the junction type (0=PARALLEL, 1=SERIES, 2=RIGID)


type = T(compIndex).type;

if type == 0 %PARALLEL
    adaptValue=0;
    for i=1:numCompEdges
        if T(compIndex).edges(i) ~= lastParentEdge %get all the others
            element_mIndex2 = T(compIndex).edges(i)+1;
            adaptValue=adaptValue+(1/Z(element_mIndex2));
        end
    end
    Z(element_mIndex)= 1/adaptValue;

    %Generate scattering matrix
    Z_values = Z(T(compIndex).edges+1); %Get the values for current comp. Z
    G_values = 1./Z_values;

    S=(2/sum(G_values))*ones(numCompEdges, 1)*G_values'-eye(numCompEdges);


    S = [S, zeros(numCompEdges, size(overallS,2)-numCompEdges)]; %zero-pad
    overallS(T(compIndex).edges+1, :) = S;

    %Save the scattering row for the parent edge, going up
    T(compIndex).scatteringUp=overallS(element_mIndex, :);

elseif type == 1 %SERIES
    adaptValue=0;
    for i=1:numCompEdges
        if T(compIndex).edges(i) ~= lastParentEdge %get all the others
            element_mIndex2 = T(compIndex).edges(i)+1;
            adaptValue=adaptValue+Z(element_mIndex2);
        end
    end
    Z(element_mIndex)= adaptValue;

    Z_values = Z(T(compIndex).edges+1); %Get the values for current comp. Z

    S=eye(numCompEdges)-(2/sum(Z_values))*Z_values(:)*ones(1, numCompEdges);


    S = [S, zeros(numCompEdges, size(overallS,2)-numCompEdges)]; %zero-pad
    overallS(T(compIndex).edges+1, :) = S;

    %Save the scattering row for the parent edge, going up
    T(compIndex).scatteringUp=overallS(element_mIndex, :);

elseif type ==2 %RIGID

    %% Rigid component adaptation



    edges = T(compIndex).edges;
    parentEdgeEndpoints = endpoints(element_mIndex, :)+1;
    edges = edges(edges ~= lastParentEdge); %remove parentedge


    %Create graph object

    endpointsLocal = endpoints(edges+1, :)+1; %offset also indexes of nodes (as needed by matlab)
    allNodes = unique(endpointsLocal); %all nodes list
    [allNodes, nodesI] = sort(allNodes);
    ids = edges+1;
    %It's important to ALWAYS pass the endpoints as STRINGS and not
    %integers when creating a graph. If we use integers MATLAB interpret
    %them as node indexes and not node names, and it will "fill" the
    %missing indexes creating spurious unconnected nodes in the graph
    edgetable = table(string(endpointsLocal), ids', 'VariableNames',{'EndNodes', 'Id'});
    %nodeTable = table(string(allNodes),'VariableNames',{'Name'});
    componentGraph = graph(edgetable);

    thevImp = adaptRJunction(componentGraph, Z(componentGraph.Edges.Id), parentEdgeEndpoints, allNodes);

    Z(element_mIndex)=thevImp;


    %% Compute scattering matrix for R node

    edges = T(compIndex).edges+1;

    %Add the virtual edge to the previous graph, since now we need
    %the complete version
    edgetable = table(string(parentEdgeEndpoints), lastParentEdge+1, 'VariableNames',{'EndNodes', 'Id'});
    componentGraph = addedge(componentGraph, edgetable);

    [B, Q, orderedEdges] = getBQ(componentGraph);

    Z_values = Z(orderedEdges);
    Z_m = diag(Z_values);

    S = getScatteringMatrix(B, Q, Z_m);

    positions = zeros(numCompEdges, 1);
    for h=1:numCompEdges
        positions(h) = find(orderedEdges == edges(h));
    end
    S=S(positions, positions); %re-order S

    S = [S, zeros(numCompEdges, size(overallS,2)-numCompEdges)]; %zero-pad
    overallS(T(compIndex).edges+1, :) = S;

    %Save the scattering row for the parent edge, going up
    T(compIndex).scatteringUp=overallS(element_mIndex, :);

       

end






end