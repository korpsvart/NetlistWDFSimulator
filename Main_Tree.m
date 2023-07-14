close all
clearvars 
clc
addpath utils\
addpath triconnectedDecomp\

%% Variables to adjust depending on the netlist

netlistFilename = 'BridgeTSP';
refEdgeId = "Vin"; %id of the edge corresponding to non-adaptable element
%Specify the names of the ports to compute the output
outputPortsIds = ["R2"];
%Specify the reference signals filenames for results validation, if you
%have any (for example output from LTSpice). Leave empty if not used
referenceSignalFilenames = ["data/audio/bridget_vr2p.wav"];
numOutputs = numel(outputPortsIds);



%% Import Input Audio Signal
[Vin,Fs] = audioread('data/audio/ExpSweep.wav');
Ts=1/Fs;

[Tree,E, Z, S, numComps, outputPorts] = parseWDFTree(netlistFilename, refEdgeId, numOutputs, outputPortsIds, Fs);


%% Simulate

[VOut] = simulateWDFTree(Tree, E, Z, S, Vin, numOutputs, numComps, outputPorts);


%% Plotting the results

numSamples = numel(Vin);


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
    xlabel('time [seconds]','Fontsize',22,'interpreter','latex');
    ylabel('$V_{\mathrm{OutR_2}}$ [V]','Fontsize',22,'interpreter','latex');
    xlim([0 tWdf(end)]);
    legend(legends, "Fontsize",22,"interpreter","latex");
    title('Output Signals','Fontsize',24,'interpreter','latex');
end

%Plot error signal

plot(tWdf,VOut(1, :)-referenceSignal(1, :),'k','Linewidth',1);grid on; xlim([0,tWdf(end)]);
xlabel('time [seconds]','Fontsize',22,'interpreter','latex');
ylabel('$E_{\mathrm{OutR_2}}$ [V]','Fontsize',22,'interpreter','latex');
title(['Error Signal'],'Fontsize',18,'interpreter','latex');


%%Compute RMSE
rmse(VOut(1, :), referenceSignal(1, :))








