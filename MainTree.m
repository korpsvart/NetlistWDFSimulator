close all
clearvars 
clc
addpath utils\common\
addpath utils\mainTree\
addpath triconnectedDecomp\

%% Variables to adjust depending on the netlist

makeGeneratorsReal = true;
netlistFilename = 'Diodo';
%Id of the edge corresponding to non-adaptable element
%(You MUST specify a non-adaptable element in this version at the moment)
%DO NOT specify a voltage generator if makeGeneratorsReal is true!
%(Output will be wrong)
refEdgeId = "D1";
%Specify the names of the ports to compute the output
outputPortsIds = ["R1"];
%Specify the reference signals filenames for results validation, if you
%have any (for example output from LTSpice). Leave empty if not used
referenceSignalFilenames = ["data/audio/diodo_output_r1_ideal.wav"];
numOutputs = numel(outputPortsIds);



%% Import Input Audio Signal
[Vin,Fs] = audioread('data/audio/ExpSweep.wav');
Ts=1/Fs;

[Tree,E, Z, S, numComps, outputPorts] = parseWDFTree(netlistFilename, refEdgeId, numOutputs, outputPortsIds, Fs, makeGeneratorsReal);


%% Simulate

[VOut] = simulateWDFTree(Tree, E, Z, S, Vin, numOutputs, numComps, outputPorts);


%% Plotting the results


plotResults(VOut, referenceSignalFilenames, Fs);









