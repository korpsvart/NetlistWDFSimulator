function [a] = get_a_element(element, b,sampleIndex, inputSignal)

type = element.Type;
switch type
    case 'R'    
        a = 0;
    case 'C'    
        a = b;
    case 'L'
        a = -b;
    case 'V'
        a = inputSignal(sampleIndex); %small series resistance value
    case 'I'
        a = 10e9*inputSignal(sampleIndex); % large resistance value
end

end

