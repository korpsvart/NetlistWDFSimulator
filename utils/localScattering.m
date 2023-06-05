function [a] = localScattering(element, b,sampleIndex, inputSignal)


type = element.Type;
switch type
    case 'V'
        a = 2*inputSignal(n)-b; %Ideal voltage generator
end

end