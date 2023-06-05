function [a] = get_a_element(element, b,sampleIndex, inputSignal)

type = element.Type;
switch type
    case 'R'    
        a = 0;
    case 'C'    
        a = b;
    case 'L'
        a = -b;
    case 'Vreal'
        a = inputSignal(sampleIndex); %small series resistance value
    case 'Ireal'
        a = 10e9*inputSignal(sampleIndex); % large resistance value
    case 'V'
        a = 2*inputSignal(sampleIndex)-b; %Ideal voltage generator
end

end

