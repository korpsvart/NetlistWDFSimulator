function [Z, S] = getZS(netlistFilename,B,Q, orderedEdges, Fs_signal, G, refEdgeId)

parsingResult = strcat(netlistFilename, '_Z_S.mat');

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
    refEdgeEndpoints = convertCharsToStrings(endpoints(refEdgeIndex, :));

    Gprime = rmedge(G, refEdgeIndex);


    graphNodesOrder = convertCharsToStrings(Gprime.Nodes.Variables);

    %Reorder Z_adapted according to graph edges ordering
    Z_reordered = zeros(size(Z_adapted, 1), 1);
    for h=1:size(Z_adapted, 1)
        Z_reordered(h)= Z_adapted(Gprime.Edges.Id(h)==adaptableEdgesIds);
    end

    
    [thevImp] = adaptRJunction(Gprime, Z_reordered, refEdgeEndpoints, graphNodesOrder);

    Z_diag(~adaptableEdgesIndexes)=thevImp;
    Z = diag(Z_diag);


    %%%%%%%%%%%%%%%%%%%

    S = getScatteringMatrix(B, Q, Z);
    
    
    Fs = Fs_signal;
    save(parsingResult,'S','Z', 'Fs');
end

    
end
    

