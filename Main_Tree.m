close all
clearvars 
clc
addpath utils\
addpath triconnectedDecomp\

%% Variables to adjust depending on the netlist


makeGeneratorsReal = false;
netlistFilename = 'BridgeTSP';
%Id of the edge corresponding to non-adaptable element
%(You MUST specify a non-adaptable element in this version at the moment)
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

[Tree,E, Z, S, numComps, outputPorts] = parseWDFTree(netlistFilename, refEdgeId, numOutputs, outputPortsIds, Fs);


%% Simulate

[VOut] = simulateWDFTree(Tree, E, Z, S, Vin, numOutputs, numComps, outputPorts);


%% Plotting the results


plotResults(VOut, referenceSignalFilenames, Fs);









