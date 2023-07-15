close all
clearvars 
clc
addpath utils\common\
addpath utils\main\


%% Variables to adjust depending on the netlist


makeGeneratorsReal = false;

netlistFilename = 'BridgeTSP';
%Id of the edge corresponding to non-adaptable element
%leave empty string if all elements are adaptable
%(including ideal generators if makeGeneratorsReal is set to true)
%DO NOT specify a voltage generator if makeGeneratorsReal is true!
%(Output will be wrong)
refEdgeId = "Vin"; 
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

[orderedEdges, Z, S] = parseWDF(netlistFilename, Fs, refEdgeId, makeGeneratorsReal);


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








 
 
