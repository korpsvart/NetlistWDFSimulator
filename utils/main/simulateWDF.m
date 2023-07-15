function [V,I] = simulateWDF(orderedEdges, Z, S, refEdgeIndex, Vin)

%% Loop a,b,I/O declaration and initialization

n = size(orderedEdges, 1);

Nsamp=length(Vin);
dimS = size(S);

ii=1;

a   = zeros(dimS(2), 1);
b   = zeros(dimS(1), 1);

%Outputs

V = zeros(n, Nsamp);
I = zeros(n, Nsamp);

%% Real-time filtering using WDF

Z_diag = diag(Z, 0);

types = orderedEdges(:, 1);

funcs = getElementsFunctions(types);

k=1:n;


if isempty(refEdgeIndex)
    %Simulate structure where all elements are adapted

     while (ii<Nsamp)
    
        %Waves reflected from adaptable elements
    
        for i=k
            a(i) = funcs{i}(b(i), Vin, ii);
        end
       
        b = S*a; %reflecting from the junction
        
        V(:, ii) = (a+b)/2;
        I(:, ii) = (a-b)./(2*Z_diag);
        
        ii = ii+1;  
     end

else
    %Simulate a structure with one non-adaptable element
    %(forward - local - backward scan)


    k=k(k~=refEdgeIndex);
    while (ii<Nsamp)
    
        %forward scan
    
        for i=k
            a(i) = funcs{i}(b(i), Vin, ii);
        end
       
    
        b = S*a; %reflecting from the junction
    
        %local scattering
        a(refEdgeIndex) = funcs{refEdgeIndex}(b(refEdgeIndex), Vin, ii);
    
        %backward scan
        b = S*a;
        
        %compute output voltages and currents
        
        V(:, ii) = (a+b)/2;
        I(:, ii) = (a-b)./(2*Z_diag);
        
        ii = ii+1;  
    end

end






end

