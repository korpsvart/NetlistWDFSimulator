close all
clearvars 
clc
addpath utils

tic

%% Variables to adjust depending on the netlist

netlistFilename = 'BridgeTSP';
refEdgeId = "Vin"; %id of the edge corresponding to non-adaptable element
%Specify the names of the ports to compute the output
outputPorts = ["R2", "R4"];
%Specify the reference signals filenames for results validation, if you
%have any (for example output from LTSpice). Leave empty if not used
referenceSignalFilenames = ["data/audio/bridget_vr2p.wav", "data/audio/bridget_vr2p.wav"];
numOutputs = numel(outputPorts);


%% Parse topology
parsingResult = strcat(netlistFilename, '.mat');

[B,Q,G,orderedEdges] = parseTopology(netlistFilename);

%% Import Input Audio Signal
[Vin,Fs] = audioread('data/audio/ExpSweep.wav');
Ts=1/Fs;

[Z, S] = getZS(netlistFilename,B,Q,orderedEdges,Fs, G,refEdgeId);


toc



%% Loop a,b,I/O declaration and initialization

n = size(orderedEdges, 1);

Nsamp=length(Vin);
dimS = size(S);

ii=1;

a   = zeros(dimS(2), 1);
b   = zeros(dimS(1), 1);

%Outputs

V = zeros(n, Nsamp);
I = zeros(n, Nsamp);

%% Real-time filtering using WDF

Z_diag = diag(Z, 0);

types = orderedEdges(:, 1);

funcs = cell(n,1);
for i=1:n
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

tic
while (ii<Nsamp)

    %forward scan

    for i=1:n
        a(i) = funcs{i}(b(i), Vin, ii);
    end
    
    %backward scan

    b = S*a; %reflecting coefficients from the junction
    
    %compute output voltages and currents
    
    V(:, ii) = (a+b)/2;
    I(:, ii) = (a-b)./(2*Z_diag);
    
    ii = ii+1;  
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

ids = orderedEdges(:, 2);
VOut=zeros(numOutputs, Nsamp);
for i=1:numOutputs
    VOut(i, :) = V(ids==outputPorts(i), :);
end

tWdf = 1/Fs*[1:Nsamp];
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








 
 
