function [] = plotResults(VOut, referenceSignalFilenames, Fs)

Nsamp = size(VOut, 2);
numOutputs = size(VOut, 1);

%Plot results of the simulation
if (~isempty(referenceSignalFilenames))
    for k=1:numel(referenceSignalFilenames)
        referenceSignal(k, :) = audioread(referenceSignalFilenames(k));
    end
    tReference = 1/Fs*[1:size(referenceSignal, 2)];
    legends = ["Reference","WDF"];
else
    legends = ["WDF"];
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


%%Compute RMSE (function only exists from MATLAB 2023)
%rmse(VOut(1, :), referenceSignal(1, :))


end

