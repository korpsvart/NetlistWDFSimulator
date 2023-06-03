close all
clearvars 
clc
addpath utils\
addpath TriconnectedDecomp\



%% Parsing netlist and building graph G

netlistFilename = 'BridgeT';


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

%% Triconnected components decomposition


%Prepare data for triconnected component decomposition
E = table2struct(G.Edges);
N = table2struct(G.Nodes);
for i=1:numel(E)
    E(i).EndNodes = convertCharsToStrings(E(i).EndNodes);
end
N = string(struct2cell(N));
[T, numComps, endpoints] = SQPRTree(E, N); %Call C++ interface
T = T(1:numComps);

%% Import Input Audio Signal
[Vin,Fs] = audioread('data/audio/ExpSweep.wav');
Ts=1/Fs;

%% Adaptation phase
numEdges = size(E, 1);
refEdgeIndex = 5; %You should give the edge of the non-linear element

[Tree, Z] = buildTree(T, numEdges, refEdgeIndex, E, Fs, endpoints);


%% Simulation phase

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
            b(edge+1)=2*Vin(n)-b(edge+1);
            a(edge+1)=b(edge+1);
        end
    end

    %Backward scan (growing depth)
    for i=1:numComps
        component = Tree(i);
        edges = component.edges;
        in = a(edges+1);
        %[sortedEdges, I] = sort(edges);
        for j=1:numel(edges)
            edge = edges(j);
            if edge~=component.parentEdge
                b(edge+1)=component.scattering(j, :)*in(:);
                if edge>=numEdges %assign only if edge is virtual
                    %Notice this check is needed if we want the output
                    %reading to be correct
                    a(edge+1)=b(edge+1);
                end
            end
        end
    end

    %Read output
    VLow(n)=(a(5)+b(5))/2;
end

%% Plotting the results

spiceOutLow = audioread('data/audio/bridget_vr2.wav');

tSpice = 1/Fs*[1:length(spiceOutLow)];
tWdf = 1/Fs*[1:numSamples];
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(tSpice,spiceOutLow,'r','Linewidth',2); hold on;
plot(tWdf, VLow,'b--','Linewidth',1); grid on;  
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
xlim([0 tSpice(end)]);
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
% plot(tSpice,spiceOutMid,'r','Linewidth',2); hold on;
% plot(tWdf, VMid,'b--','Linewidth',1); grid on; 
% xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
% ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
% xlim([0 tSpice(end)]);
% subplot(313)
% plot(tSpice,spiceOutHigh,'r','Linewidth',2); hold on;
% plot(tWdf, VHigh,'b--','Linewidth',1); grid on; 
% xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
% ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');
% xlim([0 tSpice(end)]);








