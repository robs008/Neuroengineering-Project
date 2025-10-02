function plot_parameters(data,n_canali)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot a specific parameter of the sample
% Input:        data --> Data struct
%               n_channels --> Vector containing the indices of the channels to be
%                               displayed
% Authors:  Roberto Pilotto
%           Salvatore Rapisarda
%           Chiara Razzuoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
sub=1
amp_odd=data(sub).v_eeg.V_EEG.odd_maxj(n_canali);
amp_std=data(sub).v_eeg.V_EEG.std_maxj(n_canali);

lat_odd=data(sub).v_eeg.V_EEG.odd_mean_latency(n_canali);
lat_std=data(sub).v_eeg.V_EEG.std_mean_latency(n_canali);
for sub=2:40
    amp_odd=[amp_odd,data(sub).v_eeg.V_EEG.odd_maxj(n_canali)];
    amp_std=[amp_std,data(sub).v_eeg.V_EEG.std_maxj(n_canali)];
    
    lat_odd=[lat_odd,data(sub).v_eeg.V_EEG.odd_mean_latency(n_canali)];
    lat_std=[lat_std,data(sub).v_eeg.V_EEG.std_mean_latency(n_canali)];
end

amp_odd_m=mean(amp_odd')';
amp_std_m=mean(amp_std')';

lat_odd_m=mean(lat_odd')';
lat_std_m=mean(lat_std')';

% errors
amp_odd_err=std(amp_odd')'/sqrt(40);
amp_std_err=std(amp_std')'/sqrt(40);

lat_odd_err=std(lat_odd')'/sqrt(40);
lat_std_err=std(lat_std')'/sqrt(40);


% plot
x_amp=[amp_odd_m,amp_std_m];
error_amp=[amp_odd_err,amp_std_err];


x_lat=[lat_odd_m,lat_std_m];
error_lat=[lat_odd_err,lat_std_err];



xTickLabels=convertCharsToStrings(data(1).v_eeg.V_EEG.label{n_canali(1),1});
if length(n_canali)>1
    for k=2:length(n_canali)
        xTickLabels=[xTickLabels;convertCharsToStrings(data(1).v_eeg.V_EEG.label{n_canali(k),1})];
    end
end
figure ()
subplot (2,1,1)
    b_1=bar(x_amp,'grouped');
    hold on;
    [ngroups,nbars] = size(x_amp);
    y = nan(nbars, ngroups);
    for i = 1:nbars
        y(i,:) = b_1(i).XEndPoints;
    end
    errorbar(y',x_amp,error_amp,'k','linestyle','none');
    ylabel('Amplitude (uV)')
    title('Mean Amplitude '); legend('odd','std'); 
    set(gca, 'XTick', 1:length(xTickLabels), 'XTickLabel', xTickLabels);
    
subplot (2,1,2)
    b_2=bar(x_lat,'grouped');
    hold on;
    [ngroups,nbars] = size(x_lat);
    y = nan(nbars, ngroups);
    for i = 1:nbars
        y(i,:) = b_2(i).XEndPoints;
    end
    errorbar(y',x_lat,error_lat,'k','linestyle','none');
    ylabel('Latency (s)')
    title('Mean latency '); legend('odd','std');
    set(gca, 'XTick', 1:length(xTickLabels), 'XTickLabel', xTickLabels);

end


