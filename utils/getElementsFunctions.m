function [funcs] = getElementsFunctions(types)

% Return an array of function handles
% describing the scattering relation of the elements,
% according to the types contained in the input vector

n = size(types, 1);

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

end

