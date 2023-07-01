close all
clearvars 
clc
addpath utils\
addpath triconnectedDecomp\

%% Variables to adjust depending on the netlist

netlistFilename = 'BridgeTSP';
%refEdgeIndex is the index of the edge corresponding to non-adaptable element
%(Starting from 1 in the netlist)
refEdgeId = "Vin";
%Specify the indexes of the ports you want to compute the output for
%(Starting from 1 in the netlist)
outputPortsIds = ["R2", "R4"];
%Specify the reference signals filenames for results validation, if you
%have any (for example output from LTSpice). Leave empty if not used
referenceSignalFilenames = ["data/audio/bridget_vr2p.wav", "data/audio/bridget_vr2p.wav"];
numOutputs = numel(outputPortsIds);




%% Parsing netlist and building graph G


% Parsing from LTSpice netlist

%CommentStyle option to ignore SPICE directives
%Range = 2 to ignore first line (specifies Netlist path)

M = readmatrix(strcat('data/netlist/', netlistFilename , '.txt'), 'OutputType', 'string',...
    'CommentStyle', {'.'}, 'Range', 2);

% Creating circuit graph
endNodes = M(:, 2:3);
types = extractBetween(M(:, 1), 1, 1);
ids = M(:, 1);
values = M(:, 4);
EdgeTable = table(endNodes, types, ids, values,...
'VariableNames',{'EndNodes', 'Type', 'Id', 'Value'});
G = graph(EdgeTable);

refEdgeIndex = find(EdgeTable.Id==refEdgeId);
outputPorts = zeros(1, numOutputs);
for i=1:numOutputs
    outputPorts(i) = find(EdgeTable.Id==outputPortsIds(i));
end


%% Triconnected components decomposition


%Prepare data for triconnected component decomposition
E = table2struct(EdgeTable);
N = table2struct(G.Nodes);
for i=1:numel(E)
    E(i).EndNodes = convertCharsToStrings(E(i).EndNodes);
end
N = string(struct2cell(N));
[T, numComps, endpoints] = TriconnectedComponents(E, N); %Call C++ interface
T = T(1:numComps);

%% Import Input Audio Signal
[Vin,Fs] = audioread('data/audio/ExpSweep.wav');
Ts=1/Fs;

%% Adaptation phase
numEdges = size(E, 1);

tic
[Tree, Z, S] = getZSFromTriconnected(T, numEdges, refEdgeIndex, E, Fs, endpoints);
toc

%% Simulation phase

%sort by increasing depth
[B,I] = sort([Tree.depth]);
Tree = Tree(I);

%Allocate incident and reflected waves vectors
M = size(Z, 1);
b = zeros(M, 1);
a = zeros(M, 1);


%% Setup scattering relations for real elements

types = [E.Type]';

funcs = cell(numEdges,1);
for i=1:numEdges
switch types(i)
    case 'R'    
        funcs{i} = @(b, Vin, ii) 0;
    case 'C'    
        funcs{i} = @(b, Vin, ii) b;
    case 'L'
        funcs{i} = @(b, Vin, ii) -b;
    case 'Vreal'
        funcs{i} = @(b, Vin, ii) Vin(ii); %small series resistance value
    case 'Ireal'
        funcs{i} = @(b, Vin, ii) 10e9*Vin(ii); % large resistance value
    case 'V'
        funcs{i} = @(b, Vin, ii) 2*Vin(ii)-b; % ideal voltage source
end
end




numSamples = numel(Vin);
VOut = zeros(numOutputs, numSamples);
maxDepth = Tree(end).depth;


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

tic

for n=1:numSamples
    %Forward scan
    for i=numComps:-1:1
        component = Tree(i);
        edges = component.edges;
        %Compute elements reflected waves
        %(manage linear elements)
        %The second check (with parentEdge) is only needed for the
        %particular case where edge correspond to a non-adaptable element
        %realEdges = edges(edges<numEdges & edges~=component.parentEdge);
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
toc

%% Plotting the results


if (~isempty(referenceSignalFilenames))
    for k=1:numel(referenceSignalFilenames)
        referenceSignal(k, :) = audioread(referenceSignalFilenames(k));
    end
    tReference = 1/Fs*[1:size(referenceSignal, 2)];
    legends = ["Reference","WDF"];
else
    legends = ["WDF"];
end


tWdf = 1/Fs*[1:numSamples];
figure
set(gcf, 'Color', 'w');
for i=1:numOutputs
    subplot(numOutputs, 1, i)
    if (~isempty(referenceSignalFilenames))
        plot(tReference,referenceSignal(i, :),'r','Linewidth',2); hold on;
    end
    plot(tWdf, VOut(i, :),'b--','Linewidth',1); grid on;
    xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
    ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
    xlim([0 tWdf(end)]);
    legend(legends, "Fontsize",16,"interpreter","latex");
    title('Output Signals','Fontsize',18,'interpreter','latex');
end









