function [Z] = getZElement(element, Fs)

%Get adaptation for one single element

type = element.Type;
valueString = element.Value;
value = eng2num(convertStringsToChars(valueString));
if isnan(value)
    %Quick dirty fix for unsupported input formats
    %(e.g. when netlist specify a voltage signal)
    value = 0;
end
switch type
    case 'R'    
        Z = value;
    case 'C'    
        Z = 1/(Fs*2*value);
    case 'L'
        Z = Fs*2*value;
    case 'Vreal'
        Z = 10e-6; % small series resistance value
    case 'Ireal'
        Z = 10e9; % large resistance value
end


end

