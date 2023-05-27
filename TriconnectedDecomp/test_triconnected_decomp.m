%%Code for testing the triconnected components decomposition interface




E = table2struct(G.Edges);
N = table2struct(G.Nodes);
for i=1:numel(E)
    E(i).EndNodes = convertCharsToStrings(E(i).EndNodes);
end
N = string(struct2cell(N));
[T, numComps] = SQPRTree(E, N);
T = T(1:numComps);
numEdges = size(E, 1);
refEdgeIndex = 0; %chosen arbitrarily here, in practice you should give the edge of the non-linear element

%Also pass the element values for adaptation


%NOTE: we need to add a variable parentEdge, which stores the edge
%used to reach that element. I think this is necessary for the
%simulation part
[Tree, Z] = buildTree(T, numEdges, refEdgeIndex, E, Fs);

%sort by increasing depth
[B,I] = sort([Tree.depth]);
Tree = Tree(I);

%Allocate incident and reflected waves vectors
M = size(Z, 1);
b = zeros(M, 1);
a = zeros(M, 1);


numSamples = 100;
for n=1:numSamples
    %Forward scan
    for i=numComps:-1:1
        component = Tree(i);
        edges = component.edges;
        %Compute elements reflected waves
        %(manage linear elements)
        for j=1:numel(edges)
           edge = edges(j);
           if edge < numEdges %edge is real
               element = E(edge+1);
               a(edge+1) = get_a_element(element, b(edge+1), n, Vin);
           end
        end
        %Compute junction reflected waves
        %(the actual forward scan)
        for j=1:numel(edges)
           edge = T(compIndex).edges(j);
           if edge >= numEdges %edge is virtual
               %Take all other incoming waves
               in = a(edges(edges~=edge)+1); %kinda hard to read

               %For the moment we store an individual scattering matrix
               %for each component, and the row associated is the one
               %given by the ordering of the edges inside the component.
               %It would be nice to have a single "giant" matrix S
               %where each row corresponds to an edge, but this would
               %require all rows to have same length (which is possible
               %only if we find a way to make the algorithm return only
               %3 elements components)
               b(edge+1)=component.scattering(j, :)*in;
               a(edge+1)=b(edge+1); %probably could be avoided, but it's easier to follow

           end
        end

        %Root scattering
        if i==1
            %for the time being we assume there's an ideal voltage gen
            b(refEdgeIndex+1)=2*Vin(n)-a(edge+1);

        end
    end

    %Backward scan (growing depth)
    for i=numComps:1:numComps

    end

end








