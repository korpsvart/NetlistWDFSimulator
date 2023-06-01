function [T, Z] = exploreComponent(T, N, numEdges, compIndex, lastVirtualEdge, E, Z, Fs, depth, endpoints)

%lastVirtualEdge is the virtual edge from which we came to this
%component (for the root it's the reference edge, so in theory the
%non-linear element)

T(compIndex).depth = depth;
T(compIndex).parentEdge = lastVirtualEdge;

%%Exploration phase
numCompEdges=numel(T(compIndex).edges);
for i=1:numCompEdges
    
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
for i=1:numCompEdges
    %Check if edge is real and if it's not the root non-adaptable element
    if T(compIndex).edges(i) < numEdges && T(compIndex).edges(i)~=lastVirtualEdge
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
    adaptValue=0;
    for i=1:numCompEdges
        if T(compIndex).edges(i) ~= lastVirtualEdge %get all the others
            element_mIndex2 = T(compIndex).edges(i)+1;
            adaptValue=adaptValue+(1/Z(element_mIndex2));
        end
    end
    Z(element_mIndex)= 1/adaptValue;

    %Generate scattering matrix
    Z_values = Z(T(compIndex).edges+1); %Get the values for current comp. Z
    G_values = 1./Z_values;

    S=(2/sum(G_values))*ones(numCompEdges, 1)*G_values'-eye(numCompEdges);

    %Alternatively
    Z_m = diag(Z_values);
    Q = ones(1, numCompEdges);
    %S=2*Q'*pinv(Q*inv(Z_m)*Q')*Q*inv(Z_m)-eye(numCompEdges);


    T(compIndex).scattering=S;
elseif type == 1 %SERIES
    adaptValue=0;
    for i=1:numCompEdges
        if T(compIndex).edges(i) ~= lastVirtualEdge %get all the others
            element_mIndex2 = T(compIndex).edges(i)+1;
            adaptValue=adaptValue+Z(element_mIndex2);
        end
    end
    Z(element_mIndex)= adaptValue;

    Z_values = Z(T(compIndex).edges+1); %Get the values for current comp. Z

    S=eye(numCompEdges)-(2/sum(Z_values))*Z_values(:)*ones(1, numCompEdges);

    %Alternatively

    Z_m = diag(Z_values);
    B = ones(1, numCompEdges);
    %S=eye(numCompEdges)-2*Z_m*B'*pinv(B*Z_m*B')*B;


    T(compIndex).scattering=S;
elseif type ==2 %RIGID

    %% Rigid component adaptation



    edges = T(compIndex).edges;
    edges = edges(edges < numEdges & edges ~= lastVirtualEdge); %remove parentedge


    %Create graph object

    endpointsLocal = endpoints(edges+1, :)+1; %offset also indexes of nodes (as needed by matlab)
    ids = [E.Id];
    ids = ids(endpointsLocal);
    edgetable = table(endpointsLocal, ids, 'VariableNames',{'EndNodes', 'Id'});
    componentGraph = graph(edgetable);





    %Remove one of the endpoints of the virtual edge
    endpointToRemove = endpoints(element_mIndex, 1)+1; %Take the first endpoint, for example
    referenceEndpoint = endpoints(element_mIndex, 2)+1; %Other become reference

    %We can't use the rmnode function because it will also delete incident
    %edges, which is not what we want
    %Instead we must remove the corresponding row. Hopefully matlab sorts
    %the row according to node indexes, so we can exploit the natural index
    %ordering...



    %Get adjacency matrix
    A = full(incidence(componentGraph));
    A(endpointToRemove, :) = []; %remove row


    %Get G conductance vector (and diag matrix)
    %We've the Z vector already at our disposal so it's easy
    G_vector = 1./Z(edges+1);
    G_m = diag(G_vector);



    %Compute Y
    Y = A*G_m*A'; %Admittance matrix
    Imp = inv(Y); %Impedances matrix
    %Adaptation
    Z(element_mIndex)=Imp(referenceEndpoint);


    %% Compute scattering matrix for R node

    %Add the virtual edge to the previous graph, since now we need
    %the complete version
    componentGraph = addedge(componentGraph, endpoints(lastVirtualEdge+1, 1)+1, endpoints(lastVirtualEdge+1, 2)+1);

    tree = minspantree(componentGraph);
    At = full(incidence(tree));
    q = size(At, 2); % q is the number of independent voltages

    %Extract cotree

    %Find edges in the original graph which are part of cotree
    idx = ismember(componentGraph.Edges.Id, tree.Edges.Id);

    %Convert logical array into sequential indexes
    idx_n = find(idx);

    %remove edges from original
    cotree = rmedge(componentGraph, idx_n);
    Ac = full(incidence(cotree));
    p = size(Ac, 2); % p is the number of independent currents

    % Computing F
    F = pinv(At)*Ac; %F has size q x p 

    % Computing B and Q
    B = [eye(p) -F'];
    Q = [F eye(q)];

    n = size(orderedEdges, 1);
    q = size(Q, 1);
    p = size(B, 1);

    Z_values = Z(T(compIndex).edges+1); %Get the values for current comp. Z
    Z_m = diag(Z_values);

    if (q <= p)
       Z_inv = inv(Z_m);
       S = 2*Q'*inv(Q*Z_inv*Q')*Q*Z_inv - eye(n);

       %According to MATLAb it's faster and more accurate like this
       %S = 2*Q'*(Q*(Z\Q')\Q)/Z - eye(n);
     else
       S = eye(n) - 2*Z_m*B'*inv(B*Z_m*B')*B;

       %Same as above
       %S = eye(n) - 2*Z*B'*((B*Z*B')\B);

    end

    T(compIndex).scattering=S;

    

end






end