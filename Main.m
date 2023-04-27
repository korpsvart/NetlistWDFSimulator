close all
clearvars 
clc
addpath utils

tic

netlistFilename = 'Crossover_test1';
parsingResult = strcat(netlistFilename, '.mat');

[B,Q,G,orderedEdges] = parseTopology(netlistFilename);

%% Import Input Audio Signal
[Vin,Fs] = audioread('data/audio/ExpSweep.wav');
Ts=1/Fs;

[Z, S] = get_Z_S(netlistFilename,B,Q,orderedEdges,Fs);


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
    case 'V'
        funcs{i} = @(b, Vin, ii) Vin(ii); %small series resistance value
    case 'I'
        funcs{i} = @(b, Vin, ii) 10e9*Vin(ii); % large resistance value
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

spiceOutLow = audioread('data/audio/outputlow.wav');
spiceOutMid = audioread('data/audio/outputmid.wav');
spiceOutHigh = audioread('data/audio/outputhigh.wav');

ids = orderedEdges(:, 2);
VLow = -V(ids=='Rspk3', :);
VMid = V(ids=='Rspk2', :);
VHigh = V(ids=='Rspk1', :);

tSpice = 1/Fs*[1:length(spiceOutLow)];
tWdf = 1/Fs*[1:Nsamp];
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
plot(tSpice,spiceOutMid,'r','Linewidth',2); hold on;
plot(tWdf, VMid,'b--','Linewidth',1); grid on; 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
xlim([0 tSpice(end)]);
subplot(313)
plot(tSpice,spiceOutHigh,'r','Linewidth',2); hold on;
plot(tWdf, VHigh,'b--','Linewidth',1); grid on; 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');
xlim([0 tSpice(end)]);






 
 
