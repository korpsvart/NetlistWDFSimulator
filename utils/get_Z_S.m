function [Z, S] = get_Z_S(netlistFilename,B,Q, orderedEdges, Fs_signal)

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
    
    Z = getZ(orderedEdges, Fs_signal);
   
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
    

