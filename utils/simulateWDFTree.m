function [VOut] = simulateWDFTree(Tree, E, Z, S, Vin, numOutputs, numComps, outputPorts)
%% Simulation phase

%sort by increasing depth
[~,I] = sort([Tree.depth]);
Tree = Tree(I);

%Allocate incident and reflected waves vectors
M = size(Z, 1);
b = zeros(M, 1);
a = zeros(M, 1);


%% Setup scattering relations for real elements

types = [E.Type]';
numEdges = size(E, 1);

funcs = getElementsFunctions(types);

numSamples = numel(Vin);
VOut = zeros(numOutputs, numSamples);


%Precompute real and virtual edges for each component
for i=1:numComps
    component = Tree(i);
    edges = component.edges;
    edges = edges(edges~=component.parentEdge);
    indexes = edges<numEdges;
    Tree(i).realEdges = edges(indexes); %Real edges, no parent
    Tree(i).virtualEdges = edges(~indexes); %Virtual edges, no parent
    Tree(i).nonParentEdges = [Tree(i).realEdges, Tree(i).virtualEdges];
end


for n=1:numSamples
    %Forward scan
    for i=numComps:-1:1
        component = Tree(i);
        edges = component.edges;
        %Compute elements reflected waves
        %(manage linear elements)
        for j=1:numel(component.realEdges)
           edge = component.realEdges(j);
           a(edge+1) = funcs{edge+1}(b(edge+1), Vin, n);
        end
        %Compute junction reflected waves
        %(the actual forward scan)
        %The edge is the one going up (parentEdge)
        edge = component.parentEdge;
        %Take all incoming waves
        %(Easier to take also the values we don't need, anyway if the
        %scattering matrix was built correctly it will be ignored)
        in = a(edges+1);
        a(edge+1)=component.scatteringUp(1:numel(in))*in(:);

    end

    %Root scattering
    b(edge+1)=funcs{edge+1}(a(edge+1), Vin, n);
    a(edge+1)=b(edge+1);



    %Backward scan (growing depth)
    for i=1:numComps
        component = Tree(i);
        edges = component.edges;
        in = a(edges+1);
        b(component.nonParentEdges+1)=S(component.nonParentEdges+1, 1:numel(in))*in(:);
        a(component.virtualEdges+1) = b(component.virtualEdges+1);
    end

    %Read output
    VOut(:, n)=(a(outputPorts)+b(outputPorts))/2;
end




end

