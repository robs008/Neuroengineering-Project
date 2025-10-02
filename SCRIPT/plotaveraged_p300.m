function plotaveraged_p300(data,sub,fc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot the ERPs of all EEG channels for the subject sub
%   Input:  data --> Data struct
%           sub --> Subject ID (from 1 to 40)
%           fc --> Sampling frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
epoch_camp=256;
n_before=51;
figure
for i=1:30
    subplot(5,6,i)
    plot((-n_before:epoch_camp-1-n_before)/(fc),data(sub).v_eeg.V_EEG.erp_odd(i,:));
    title([data(sub).v_eeg.V_EEG.label{i,1},' ERP oddball sub ', num2str(sub)])
    xlabel('Time (s)'); ylabel('Amplitude (uV)');
end

figure
for i=1:30
    subplot(5,6,i)
    plot((-n_before:epoch_camp-1-n_before)/(fc),data(sub).v_eeg.V_EEG.erp_std(i,:));
    title([data(sub).v_eeg.V_EEG.label{i,1},' ERP standard sub ', num2str(sub)]);
    xlabel('Time (s)'); ylabel('Amplitude (uV)');
end
end