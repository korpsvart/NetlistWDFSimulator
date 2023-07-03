function [V,I] = simulate(orderedEdges, Z, S, refEdgeIndex, Vin)

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

funcs = cell(n,1);
for i=1:n
switch types(i)
    case 'R'    
        funcs{i} = @(b, Vin, ii) 0;
    case 'C'    
        funcs{i} = @(b, Vin, ii) b;
    case 'L'
        funcs{i} = @(b, Vin, ii) -b;
    case 'Vreal'
        funcs{i} = @(b, Vin, ii) Vin(ii); %small series resistance value
    case 'Ireal'
        funcs{i} = @(b, Vin, ii) 10e9*Vin(ii); % large resistance value
    case 'V'
        funcs{i} = @(b, Vin, ii) 2*Vin(ii)-b; % ideal voltage source
end
end

k=1:n;
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

