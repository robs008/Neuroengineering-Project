function time_trend(data,sub,n_canali)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot the temporal trend of P300 amplitudes and latencies for subject sub
%   Input:  data --> Data struct
%           sub --> Subject ID (from 1 to 40)
%           n_channels --> Vector containing the indices of the channels to be
%                          displayed
%   % Authors:  Roberto Pilotto
%           Salvatore Rapisarda
%           Chiara Razzuoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sub_sp=sub;
first_odd=data(1).v_eeg.V_EEG.odd_max3j(n_canali,1);
middle_odd=data(1).v_eeg.V_EEG.odd_max3j(n_canali,2);
last_odd=data(1).v_eeg.V_EEG.odd_max3j(n_canali,3);
first_std=data(1).v_eeg.V_EEG.std_max3j(n_canali,1);
middle_std=data(1).v_eeg.V_EEG.std_max3j(n_canali,2);
last_std=data(1).v_eeg.V_EEG.std_max3j(n_canali,3);
for sub=2:40
        first_odd=[first_odd,data(sub).v_eeg.V_EEG.odd_max3j(n_canali,1)];
        middle_odd=[middle_odd,data(sub).v_eeg.V_EEG.odd_max3j(n_canali,2)];
        last_odd=[last_odd,data(sub).v_eeg.V_EEG.odd_max3j(n_canali,3)];
        first_std=[first_std,data(sub).v_eeg.V_EEG.std_max3j(n_canali,1)];
        middle_std=[middle_std,data(sub).v_eeg.V_EEG.std_max3j(n_canali,2)];
        last_std=[last_std,data(sub).v_eeg.V_EEG.std_max3j(n_canali,3)];
end
first_odd_m=mean(first_odd')';
middle_odd_m=mean(middle_odd')';
last_odd_m=mean(last_odd')';
first_std_m=mean(first_std')';
middle_std_m=mean(middle_std')';
last_std_m=mean(last_std')';
clear std
first_odd_err=std(first_odd')'/sqrt(40);
middle_odd_err=std(middle_odd')'/sqrt(40);
last_odd_err=std(last_odd')'/sqrt(40);
first_std_err=std(first_std')'/sqrt(40);
middle_std_err=std(middle_std')'/sqrt(40);
last_std_err=std(last_std')'/sqrt(40);

x_odd=[first_odd_m,middle_odd_m,last_odd_m];
err_odd=[first_odd_err,middle_odd_err,last_odd_err];
x_std=[first_std_m,middle_std_m,last_std_m];
err_std=[first_std_err,middle_std_err,last_std_err];
xTickLabels=convertCharsToStrings(data(1).v_eeg.V_EEG.label{n_canali(1),1});
if length(n_canali)>1
    for k=2:length(n_canali)
        xTickLabels=[xTickLabels;convertCharsToStrings(data(1).v_eeg.V_EEG.label{n_canali(k),1})];
    end
end
figure ()
subplot (2,1,1)
    b_1=bar(x_odd,'grouped');
    hold on;
    [ngroups,nbars] = size(x_odd);
    y = nan(nbars, ngroups);
    for i = 1:nbars
        y(i,:) = b_1(i).XEndPoints;
    end
    errorbar(y',x_odd,err_odd,'k','linestyle','none');
    ylabel('Amplitude (uV)')
    title('Amplitude oddball '); legend('First stimuli','Middle stimuli','Last stimuli'); 
    set(gca, 'XTick', 1:length(xTickLabels), 'XTickLabel', xTickLabels);
    
subplot (2,1,2)
    b_2=bar(x_std,'grouped');
    hold on;
    [ngroups,nbars] = size(x_std);
    y = nan(nbars, ngroups);
    for i = 1:nbars
        y(i,:) = b_2(i).XEndPoints;
    end
    errorbar(y',x_std,err_std,'k','linestyle','none');
    ylabel('Amplitude (uV)')
    title('Amplitude standard '); legend('First stimuli','Middle stimuli','Last stimuli');
    set(gca, 'XTick', 1:length(xTickLabels), 'XTickLabel', xTickLabels);

sub=sub_sp;
x_odd=mean(data(sub).v_eeg.V_EEG.epoch_odd_latency(n_canali,:));
y_odd=std(data(sub).v_eeg.V_EEG.epoch_odd_latency(n_canali,:))/sqrt(length(n_canali));

x_std=mean(data(sub).v_eeg.V_EEG.epoch_std_latency(n_canali,:));
y_std=std(data(sub).v_eeg.V_EEG.epoch_std_latency(n_canali,:))/sqrt(length(n_canali));

figure
subplot(2,1,1)
    plot(x_odd,'b','LineWidth',2); hold on;
    plot(x_odd+y_odd,'b'); hold on;
    plot(x_odd-y_odd,'b');
    xlabel('# Event'); ylabel('Latency (s)'), title(['Time trend latency oddball subject ',num2str(sub)]);
subplot(2,1,2)
    plot(x_std,'b','LineWidth',2); hold on;
    plot(x_std+y_std,'b'); hold on;
    plot(x_std-y_std,'b');
    xlabel('# Event'); ylabel('Latency (s)'), title(['Time trend latency standard subject ',num2str(sub)]);

end
