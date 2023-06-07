function [Z] = getZ(elements, Fs)
n = size(elements, 1);
v = zeros(n, 1);

for i=1:n
    type = elements(i, 1);
    valueString = elements(i, 3);
    value = eng2num(convertStringsToChars(valueString));
    if isnan(value)
        %Quick dirty fix for unsupported input formats
        %(e.g. when netlist specify a voltage signal)
        value = 0;
    end
    switch type
        case 'R'    
            v(i) = value;
        case 'C'    
            v(i) = 1/(Fs*2*value);
        case 'L'
            v(i) = Fs*2*value;
        case 'Vreal'
            v(i) = 10e-6; % small series resistance value
        case 'Ireal'
            v(i) = 10e9; % large resistance value
    end
end

Z = diag(v);

end

