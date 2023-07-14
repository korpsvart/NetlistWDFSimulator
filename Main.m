close all
clearvars 
clc
addpath utils
addpath utils_profiling\


%% Variables to adjust depending on the netlist

netlistFilename = 'BridgeTSP';
refEdgeId = "Vin"; %id of the edge corresponding to non-adaptable element
%Specify the names of the ports to compute the output
outputPorts = ["R2"];
%Specify the reference signals filenames for results validation, if you
%have any (for example output from LTSpice). Leave empty if not used
referenceSignalFilenames = ["data/audio/bridget_vr2p.wav"];
numOutputs = numel(outputPorts);

%% Import Input Audio Signal
[Vin,Fs] = audioread('data/audio/ExpSweep.wav');
Ts=1/Fs;

%% Build WDF Structure

[orderedEdges, Z, S] = parseWDF(netlistFilename, Fs, refEdgeId);


%% Simulate

ids = orderedEdges(:, 2);
refEdgeIndex = find(ids==refEdgeId);

[V, I] = simulateWDF(orderedEdges, Z, S, refEdgeIndex, Vin);


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

Nsamp = length(Vin);
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


%%Compute RMSE
rmse(VOut(1, :), referenceSignal(1, :))








 
 
