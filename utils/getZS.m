function [Z, S] = getZS(netlistFilename,B,Q, orderedEdges, Fs_signal, G, refEdgeId)

parsingResult = strcat(netlistFilename, '_Z_S.mat');

n = size(orderedEdges, 1);
q = size(Q, 1);
p = size(B, 1);
previousParsingLoaded = false;
updated = false;

if (isfile(parsingResult))
    %Load parsing result if available
    load(parsingResult, 'Z', 'S', 'Fs');
    
    previousParsingLoaded=true;
    
    M_new = readmatrix(strcat('data/netlist/', netlistFilename , '.txt'), 'OutputType', 'string',...
    'CommentStyle', {'.'}, 'Range', 2);
    %Check for updated values and update orderedEdges if needed
    [updated, orderedEdges] = checkValueUpdates(orderedEdges, M_new);


    if (Fs_signal ~= Fs)
        %Check if an old sampling frequency was loaded
        %and if it's different from the new one
        updated = true;
    end
    
end


%% Computing Scattering Matrix
if (~previousParsingLoaded || updated)
    
    adaptableEdgesIndexes = orderedEdges(:, 2)~=refEdgeId;
    Z = getZ(orderedEdges, Fs_signal);
    adaptableEdgesIds = orderedEdges(adaptableEdgesIndexes, 2);
    Z_diag = diag(Z);
    Z_adapted = Z_diag(adaptableEdgesIndexes);


    %%%%%%%%%%%%%%%%%%
    % Compute the Z for the non-adaptable element

    %get endpoints info
    endpoints = G.Edges.EndNodes;
    refEdgeIndex = find(G.Edges.Id == refEdgeId);
    refEdgeEndpoints = endpoints(refEdgeIndex, :);
    allNodes = unique(endpoints); %all nodes list
    [allNodes, nodesI] = sort(allNodes);
    allNodes = convertCharsToStrings(allNodes);

    Gprime = rmedge(G, refEdgeIndex);

    endpointToRemove = find(allNodes == refEdgeEndpoints(1));
    referenceEndpoint =find(allNodes == refEdgeEndpoints(2));

    A = full(incidence(Gprime));
    A(endpointToRemove, :) = [];

    if (referenceEndpoint>endpointToRemove) %shift for removal
        referenceEndpoint = referenceEndpoint-1;
    end

    %Reorder Z_adapted according to graph edges ordering
    Z_reordered = zeros(size(Z_adapted, 1), 1);
    for h=1:size(Z_adapted, 1)
        Z_reordered(h)= Z_adapted(Gprime.Edges.Id(h)==adaptableEdgesIds);
    end

    G_vector = 1./Z_reordered;
    G_m = diag(G_vector);

    %Compute Y
    Y = A*G_m*A'; %Admittance matrix
    Imp = inv(Y); %Impedances matrix
    %Adaptation
    Z_diag(~adaptableEdgesIndexes)=Imp(referenceEndpoint, referenceEndpoint);
    Z = diag(Z_diag);

    %reorder columns of A
    % positions = zeros(size(A, 2), 1);
    % for h=1:size(A, 2);
    %     positions(h) = find(orderedEdges == edges(h));
    % end
    % A










    %%%%%%%%%%%%%%%%%%%


   
    if (q <= p)
       Z_inv = inv(Z);
       S = 2*Q'*inv(Q*Z_inv*Q')*Q*Z_inv - eye(n);

       %According to MATLAb it's faster and more accurate like this
       %S = 2*Q'*(Q*(Z\Q')\Q)/Z - eye(n);
     else
       S = eye(n) - 2*Z*B'*inv(B*Z*B')*B;

       %Same as above
       %S = eye(n) - 2*Z*B'*((B*Z*B')\B);

    end
    
    Fs = Fs_signal;
    save(parsingResult,'S','Z', 'Fs');
end


    
    
    
    
end
    

