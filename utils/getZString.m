function [ZString] = getZString(elements)
n = size(elements, 1);
ZString = "Z=diag([";

for i=1:n
    type = elements(i, 1);
    name = elements(i, 2);
    switch type
        case 'R'    
            ZString = append(ZString, name);
        case 'C'
            ZString = append(ZString, sprintf("Ts/(2*%s)", name));
        case 'L'
            ZString = append(ZString, sprintf("2*Fs*%s", name));
        case 'V'
            ZString = append(ZString, sprintf("%s_RSer", name)); % small series resistance value
        case 'I'
            ZString = append(ZString, sprintf("%s_RPar", name)); % large resistance value
    end
    if i < n
        ZString = append(ZString, ", ");
    end
end

ZString = append(ZString, "])");

end

