close all
clearvars 
clc
addpath utils


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

%% Build WDF Structure

[orderedEdges, Z, S] = parseWDF(netlistFilename, Fs, refEdgeId);


%% Simulate

ids = orderedEdges(:, 2);
refEdgeIndex = find(ids==refEdgeId);

[V, I] = simulateWDF(orderedEdges, Z, S, refEdgeIndex, Vin);


%% Plotting the results

Nsamp = length(Vin);
VOut=zeros(numOutputs, Nsamp);
for i=1:numOutputs
    VOut(i, :) = V(ids==outputPortsIds(i), :);
end

plotResults(VOut, referenceSignalFilenames, Fs);








 
 
