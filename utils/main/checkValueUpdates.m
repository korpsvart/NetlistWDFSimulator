function [updated, orderedEdges] = checkValueUpdates(orderedEdges, newM)
n = size(orderedEdges, 1);
updated = false;
for i=1:n
    id = orderedEdges(i, 2);
    old_value = orderedEdges(i, 3);
    new_value = newM(newM(:, 1)==id,4);
    if (old_value ~= new_value)
        updated = true; %flag as updated
        orderedEdges(i, 3) = new_value; %update
    end
end
end

