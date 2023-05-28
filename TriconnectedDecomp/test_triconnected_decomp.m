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
refEdgeIndex = 11; %chosen arbitrarily here, in practice you should give the edge of the non-linear element

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


numSamples = numel(Vin);
Vout = zeros(numSamples, 1);
maxDepth = Tree(end).depth;
for n=1:numSamples
    %Forward scan
    for i=numComps:-1:1
        component = Tree(i);
        edges = component.edges;
        %Compute elements reflected waves
        %(manage linear elements)
        %The second check (with parentEdge) is only needed for the
        %particular case where edge correspond to a non-adaptable element
        realEdges = edges(edges<numEdges & edges~=component.parentEdge);
        for j=1:numel(realEdges)
           edge = realEdges(j);
           element = E(edge+1);
           a(edge+1) = get_a_element(element, b(edge+1), n, Vin);
        end
        %Compute junction reflected waves
        %(the actual forward scan)
        %The edge is the one going up (parentEdge)
        edge = component.parentEdge;
        parentEdgeIndex = find(edges==edge);
        %Take all incoming waves
        %(Easier to take also the values we don't need, anyway if the
        %scattering matrix was built correctly it will be ignored)
        in = a(edges+1);
  
        %For the moment we store an individual scattering matrix
        %for each component, and the row associated is the one
        %given by the ordering of the edges inside the component.
        %It would be nice to have a single "giant" matrix S
        %where each row corresponds to an edge, but this would
        %require all rows to have same length (which is possible
        %only if we find a way to make the algorithm return only
        %3 elements components)
        b(edge+1)=component.scattering(parentEdgeIndex, :)*in(:);
        a(edge+1)=b(edge+1); %should be avoided for the root, but it should do no damage anyway

        %Root scattering
        if i==1
            %for the time being we assume there's an ideal voltage gen
            b(edge+1)=Vin(n)-a(edge+1);
            a(edge+1)=b(edge+1);
        end
    end

    %Backward scan (growing depth)
    for i=1:numComps
        component = Tree(i);
        edges = component.edges;
        in = a(edges+1);
        for j=1:numel(edges)
            edge = edges(j);
            if edge~=component.parentEdge
                b(edge+1)=component.scattering(j, :)*in(:);
                if component.depth<maxDepth
                    %Is this check really needed?
                    %Is it safe to assign this either way since it will
                    %simply get replaced, or does it messes up with the
                    %output reading?
                    a(edge+1)=b(edge+1);
                end
            end
        end
    end

    %Read output
    VHigh2(n)=(a(9)+b(9))/2;
    VMid2(n)=(a(6)+b(6))/2;
    VLow2(n)=(a(1)+b(1))/2;
end








