%% Neuroengineering Project-Academic Year 2023-2024.
% "EEG Signals–Active Visual Oddball Stimuli" Group 38
clc;
clear all;
close all;

%% LOADING DATA

% Initialize data struct array
data = struct('mnt', {}, 'v_eeg', {}, 'v_event', {}, 'age', [], 'gender', '', 'handedness', '');

Age = [20,24,18,21,23,22,20,25,30,20,19,25,22,30,21,20,19,23,19,20,23,20,22,21,19,20,21,21,20,19,20,21,28,18,24,29,22,19,21,21];
Gender = ['M','F','F','F','M','M','F','F','F','F','F','M','M','M','M','F','F','F','M','F','M','F','F','M','F','M','M','F','F','F','F','M','F','F','M','F','F','F','F','M'];
Handedness = ['r','r','r','r','r','r','r','r','r','r','r','l','r','r','r','r','r','r','l','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r'];

% Insert here the data folder
data_folder='C:\Users\Roberto\Desktop\Università\Neuroengineering\project\DATA';

for k = 1:40
    path = [data_folder,'\S', num2str(k), '\'];

    % Load data into struct array
    data(k).mnt = load([path, 'MNT.mat']);
    data(k).v_eeg = load([path, 'V_EEG.mat']);
    data(k).v_event = load([path, 'V_event.mat']);
    data(k).age = Age(k);
    data(k).gender = Gender(k);
    data(k).handedness = Handedness(k);
end

%% Useful data

% Sampling frequency
fc=256;

% Channels most affected by the eye blink artifact
Fp=["FP1","FP2"];
for i=1:length(Fp)
    for k=1:30
        if strcmp(Fp(i),data(1).v_eeg.V_EEG.label{k,1})==1
            Fp_chan(i)=k;
        end
    end
end


%% Beginning of struct cycle

for sub=1:40

    % HPF 0.1HZ definition
    n=4;
    wn=0.1/(fc/2);
    [b,a]=butter(n,wn,'high');
    
    % LPF 35HZ definition
    n=6;
    wn=35/(fc/2);
    [b_l,a_l]=butter(n,wn,'low');    
    
    % NaN control
    [ind_nan_r,ind_nan_c]=find(isnan(data(sub).v_eeg.V_EEG.trial{1, 1}));
    data(sub).v_eeg.V_EEG.trial{1, 1}(ind_nan_r,ind_nan_c)=0;

    % HPF
    data(sub).v_eeg.V_EEG.trial{1, 1}=filtfilt(b,a,data(sub).v_eeg.V_EEG.trial{1, 1}')';

    % % Indipendent Components Analysis ICA
    %
    % % First of all Principal Components Analysis (PCA)
    % q=21      % Number of components
    % [coeff, data_pca,latent,tsquared,explained,mu]=pca(data(sub).v_eeg.V_EEG.trial{1, 1}',NumComponents=q);
    % % ICA
    % Mdl = rica(data_pca, q)
    % Data_ICA = transform (Mdl, data_pca) ;
    %
    % In order to identify the components most affected by the eye blink
    % we are searching for the components with a higher correlation with
    % the EOG signals
    % for i=1:q
    %     cross_corr(i,1)=xcorr(Data_ICA(:,i), data(sub).v_eeg.V_EEG.trial{1, 1}(31,:),0,'normalized');
    %     cross_corr(i,2)=xcorr(Data_ICA(:,i), data(sub).v_eeg.V_EEG.trial{1, 1}(32,:),0,'normalized');
    %     cross_corr(i,3)=xcorr(Data_ICA(:,i), data(sub).v_eeg.V_EEG.trial{1, 1}(33,:),0,'normalized');
    % end
    % cross_corr=abs(cross_corr);
    % cross_corr_max=max(cross_corr')';
    % blink_components=find(cross_corr_max>0.9);
    % Data_ICA_noBlinks = Data_ICA;
    % % Removal of the identified components
    % Data_ICA_noBlinks(:,blink_components)=zeros(length(Data_ICA(:,1)),length(blink_components));
    % % Reconstruction of the signals
    % Data_PCA_noBlinks=Data_ICA_noBlinks*Mdl.TransformWeights;
    % data(sub).v_eeg.V_EEG.trial{1, 1}=(Data_PCA_noBlinks*coeff')';

    % Offline reference P9
    ref=data(sub).v_eeg.V_EEG.trial{1, 1}(9,:);
    for k=1:30
        data(sub).v_eeg.V_EEG.trial{1, 1}(k,:)=data(sub).v_eeg.V_EEG.trial{1, 1}(k,:)-ref;
    end

    % Turning off FP1 and FP2
    data(sub).v_eeg.V_EEG.trial{1, 1}(Fp_chan,:)=zeros(length(Fp_chan), length(data(sub).v_eeg.V_EEG.trial{1, 1}(1,:)));

    % LPF
    data(sub).v_eeg.V_EEG.trial{1, 1}=filtfilt(b_l,a_l,data(sub).v_eeg.V_EEG.trial{1, 1}')';

    % Initialization of onset vectors
    data(sub).v_eeg.V_EEG.onset_odd = [];
    data(sub).v_eeg.V_EEG.onset_std = [];
    data(sub).v_eeg.V_EEG.onset = [];

    % The onsets of oddball and standard events, as well as all events, are loaded
    k=1;
    a=size(data(sub).v_event.V_event);
    for i = 1:a(2)
        if strcmp(data(sub).v_event.V_event(i).type, 'trigger') && data(sub).v_event.V_event(i).value == 1 
            if data(sub).v_event.V_event(i+1).response==201
                data(sub).v_eeg.V_EEG.onset=[data(sub).v_eeg.V_EEG.onset,data(sub).v_event.V_event(i).onset];
                data(sub).v_eeg.V_EEG.gt(k)=1;
                data(sub).v_eeg.V_EEG.onset_odd = [data(sub).v_eeg.V_EEG.onset_odd, data(sub).v_event.V_event(i).onset];
                k=k+1;
            elseif data(sub).v_event.V_event(i+1).response==202
                data(sub).v_eeg.V_EEG.onset=[data(sub).v_eeg.V_EEG.onset,data(sub).v_event.V_event(i).onset];
                data(sub).v_eeg.V_EEG.gt(k)=0;
                data(sub).v_eeg.V_EEG.onset_std = [data(sub).v_eeg.V_EEG.onset_std, data(sub).v_event.V_event(i).onset];
                k=k+1;
            end
        elseif strcmp(data(sub).v_event.V_event(i).type, 'trigger') && data(sub).v_event.V_event(i).value == 2 
            if data(sub).v_event.V_event(i+1).response==201
                data(sub).v_eeg.V_EEG.onset=[data(sub).v_eeg.V_EEG.onset,data(sub).v_event.V_event(i).onset];
                data(sub).v_eeg.V_EEG.gt(k)=0;
                data(sub).v_eeg.V_EEG.onset_std = [data(sub).v_eeg.V_EEG.onset_std, data(sub).v_event.V_event(i).onset];
                k=k+1;
            elseif data(sub).v_event.V_event(i+1).response==202
                data(sub).v_eeg.V_EEG.onset=[data(sub).v_eeg.V_EEG.onset,data(sub).v_event.V_event(i).onset];
                data(sub).v_eeg.V_EEG.gt(k)=1;
                data(sub).v_eeg.V_EEG.onset_odd = [data(sub).v_eeg.V_EEG.onset_odd, data(sub).v_event.V_event(i).onset];
                k=k+1;
            end
        end
    end
    
    % 1s epoch from 200 ms before to 800 ms after the event
    epoch_camp=round(fc);
    n_before=round(epoch_camp/5);  %200ms
    n_after=epoch_camp-n_before;   %800ms

    % Search for epochs of oddball events
    for i=1:length(data(sub).v_eeg.V_EEG.onset_odd)
        in=round(data(sub).v_eeg.V_EEG.onset_odd(i)*fc);
        data(sub).v_eeg.V_EEG.epoch_odd(1:30,((i-1)*epoch_camp+1):(i*epoch_camp))=data(sub).v_eeg.V_EEG.trial{1, 1}(1:30,in-n_before:in+epoch_camp-1-n_before);
    end
    
    % Search for epochs of standard events
    for i=1:length(data(sub).v_eeg.V_EEG.onset_std)
        in=round(data(sub).v_eeg.V_EEG.onset_std(i)*fc);
        data(sub).v_eeg.V_EEG.epoch_std(1:30,((i-1)*epoch_camp+1):(i*epoch_camp))=data(sub).v_eeg.V_EEG.trial{1, 1}(1:30,in-n_before:in+epoch_camp-1-n_before);
    end

    % Search for epochs of all events
    for i=1:length(data(sub).v_eeg.V_EEG.onset)
        in=round(data(sub).v_eeg.V_EEG.onset(i)*fc);
        data(sub).v_eeg.V_EEG.epoch(1:30,((i-1)*epoch_camp+1):(i*epoch_camp))=data(sub).v_eeg.V_EEG.trial{1, 1}(1:30,in-n_before:in+epoch_camp-1-n_before);
    end
    
    % Averaging to highlight the waveforms of the signals
    for i=1:30
        data(sub).v_eeg.V_EEG.erp_odd(i,:)=av_j(data(sub).v_eeg.V_EEG.epoch_odd(i,:),length(data(sub).v_eeg.V_EEG.onset_odd),epoch_camp);
        data(sub).v_eeg.V_EEG.erp_std(i,:)=av_j(data(sub).v_eeg.V_EEG.epoch_std(i,:),length(data(sub).v_eeg.V_EEG.onset_std),epoch_camp);
        data(sub).v_eeg.V_EEG.erp(i,:)=av_j(data(sub).v_eeg.V_EEG.epoch(i,:),length(data(sub).v_eeg.V_EEG.onset),epoch_camp);
    end
    
    % Search for P300 latency
    inizio=round((0.2+0.2)*fc);
    fine=round((0.2+0.5)*fc);
    l=length(inizio:fine);
    t_ax=0:l-1;
    t_ax=t_ax/fc+0.2;
    
    for i=1:30
        for k=1:length(data(sub).v_eeg.V_EEG.epoch_odd(1,:))/epoch_camp
                [data(sub).v_eeg.V_EEG.epoch_odd_max(i,k),data(sub).v_eeg.V_EEG.epoch_odd_latency(i,k)]=max(data(sub).v_eeg.V_EEG.epoch_odd(i,(k-1)*epoch_camp+1+inizio:(k-1)*epoch_camp+1+fine));
                data(sub).v_eeg.V_EEG.epoch_odd_latency(i,k)=t_ax(data(sub).v_eeg.V_EEG.epoch_odd_latency(i,k));
        end
    end
    for i=1:30
        for k=1:length(data(sub).v_eeg.V_EEG.epoch_std(1,:))/epoch_camp
                [data(sub).v_eeg.V_EEG.epoch_std_max(i,k),data(sub).v_eeg.V_EEG.epoch_std_latency(i,k)]=max(data(sub).v_eeg.V_EEG.epoch_std(i,(k-1)*epoch_camp+1+inizio:(k-1)*epoch_camp+1+fine));
                data(sub).v_eeg.V_EEG.epoch_std_latency(i,k)=t_ax(data(sub).v_eeg.V_EEG.epoch_std_latency(i,k));
        end
    end
    for i=1:30
        for k=1:length(data(sub).v_eeg.V_EEG.epoch(1,:))/epoch_camp
                [data(sub).v_eeg.V_EEG.epoch_max(i,k),data(sub).v_eeg.V_EEG.epoch_latency(i,k)]=max(data(sub).v_eeg.V_EEG.epoch(i,(k-1)*epoch_camp+1+inizio:(k-1)*epoch_camp+1+fine));
                data(sub).v_eeg.V_EEG.epoch_latency(i,k)=t_ax(data(sub).v_eeg.V_EEG.epoch_latency(i,k));
        end
    end

    % Calculation of peaks with jitter correction, leveraging the latency estimated in the previous step
    data=max_jav(data,sub,fc,epoch_camp,n_before);
    

    data(sub).v_eeg.V_EEG.std_mean_max=mean(data(sub).v_eeg.V_EEG.epoch_std_max')';
    data(sub).v_eeg.V_EEG.std_err_max=std(data(sub).v_eeg.V_EEG.epoch_std_max')'/sqrt(length(data(sub).v_eeg.V_EEG.epoch_std_max(1,:)));
    data(sub).v_eeg.V_EEG.odd_mean_max=mean(data(sub).v_eeg.V_EEG.epoch_odd_max')';
    data(sub).v_eeg.V_EEG.odd_err_max=std(data(sub).v_eeg.V_EEG.epoch_odd_max')'/sqrt(length(data(sub).v_eeg.V_EEG.epoch_odd_max(1,:)));
    data(sub).v_eeg.V_EEG.mean_max=mean(data(sub).v_eeg.V_EEG.epoch_max')';
    data(sub).v_eeg.V_EEG.err_max=std(data(sub).v_eeg.V_EEG.epoch_max')'/sqrt(length(data(sub).v_eeg.V_EEG.epoch_max(1,:)));
    
    data(sub).v_eeg.V_EEG.std_mean_latency=mean(data(sub).v_eeg.V_EEG.epoch_std_latency')';
    data(sub).v_eeg.V_EEG.std_err_latency=std(data(sub).v_eeg.V_EEG.epoch_std_latency')'/sqrt(length(data(sub).v_eeg.V_EEG.epoch_std_latency(1,:)));
    data(sub).v_eeg.V_EEG.odd_mean_latency=mean(data(sub).v_eeg.V_EEG.epoch_odd_latency')';
    data(sub).v_eeg.V_EEG.odd_err_latency=std(data(sub).v_eeg.V_EEG.epoch_odd_latency')'/sqrt(length(data(sub).v_eeg.V_EEG.epoch_odd_latency(1,:)));
    data(sub).v_eeg.V_EEG.mean_latency=mean(data(sub).v_eeg.V_EEG.epoch_latency')';
    data(sub).v_eeg.V_EEG.err_latency=std(data(sub).v_eeg.V_EEG.epoch_latency')'/sqrt(length(data(sub).v_eeg.V_EEG.epoch_latency(1,:)));

    % % Identification of the P300 based on a threshold
    % data(sub).v_eeg.V_EEG.threshold=(data(sub).v_eeg.V_EEG.std_mean_max+data(sub).v_eeg.V_EEG.odd_mean_max)/2;
    % data(sub).v_eeg.V_EEG.prediction_th=data(sub).v_eeg.V_EEG.epoch_max > data(sub).v_eeg.V_EEG.threshold;
    % data(sub).v_eeg.V_EEG.prediction_th_cm=confusionmat(data(sub).v_eeg.V_EEG.gt,double(data(sub).v_eeg.V_EEG.prediction_th(10,:)));
    % 
    % % Identification of the P300 based on the cross-correlation coefficient
    % for i=1:30
    %     for k=1:length(data(sub).v_eeg.V_EEG.onset)
    %         %   a--> related to the average oddball ERP
    %         %   b--> related to the average standard ERP
    %         a(i,k)=xcorr( data(sub).v_eeg.V_EEG.epoch(i,(k-1)*epoch_camp+1:k*epoch_camp) , data(sub).v_eeg.V_EEG.erp_odd(i,:),0,'normalized' );
    %         b(i,k)=xcorr( data(sub).v_eeg.V_EEG.epoch(i,(k-1)*epoch_camp+1:k*epoch_camp) ,data(sub).v_eeg.V_EEG.erp_std(i,:),0,"normalized");
    %     end
    % end
    % data(sub).v_eeg.V_EEG.prediction_cc=a>b;
    % data(sub).v_eeg.V_EEG.prediction_cc_cm=confusionmat(data(sub).v_eeg.V_EEG.gt,double(data(sub).v_eeg.V_EEG.prediction_cc(10,:)));
   

    % Error analysis: check if participants make errors and, in this case, examine their response time when making errors
    % or when they do not. There might be a trend, for example: when making an error, the reaction time is shorter.
    % On the other hand, those with longer reaction times might make fewer errors.
    a=size(data(sub).v_event.V_event);
    k=1;
    for i = 1:a(2)
        if strcmp(data(sub).v_event.V_event(i).type, 'trigger') && data(sub).v_event.V_event(i).value ~= 0 
            if data(sub).v_event.V_event(i+1).response==201
                data(sub).v_eeg.V_EEG.errors(k)=0;    % no error
                k=k+1;
            elseif data(sub).v_event.V_event(i+1).response==202
                data(sub).v_eeg.V_EEG.errors(k)=1;    % error
                k=k+1;
            end
        end
    end
    for k=1:length(data(sub).v_eeg.V_EEG.errors)
        data(sub).v_eeg.V_EEG.errors_cumulative(k)=sum(data(sub).v_eeg.V_EEG.errors(1:k));
    end
end

%% Plot section

% Individual patient plots

% Channels where the P300 ERP is visually more prominent on the map
n_chan=["Pz";"PO3";"P3";"P4";"PO4"];
for i=1:length(n_chan)
    for k=1:30
        if strcmp(n_chan(i),data(1).v_eeg.V_EEG.label{k,1})==1
            n_channel(i)=k;
        end
    end
end
plot_parameters(data,n_channel)

% Statistical test
p_value=[];
sub=1;
dati_condizione_1 = data(sub).v_eeg.V_EEG.odd_maxj(n_channel,:);  
dati_condizione_2 = data(sub).v_eeg.V_EEG.std_maxj(n_channel,:);
for sub=2:40
        dati_condizione_1 = [dati_condizione_1,data(sub).v_eeg.V_EEG.odd_maxj(n_channel,:)];  
        dati_condizione_2 = [dati_condizione_2,data(sub).v_eeg.V_EEG.std_maxj(n_channel,:)]; 
end
for k=1:length(n_channel)
        [~, p_value(k)] = ttest(dati_condizione_1(k,:), dati_condizione_2(k,:));
end

time_trend(data,37,n_channel)

% Plots related to all patients

odd=data(1).v_eeg.V_EEG.odd_maxj;
std=data(1).v_eeg.V_EEG.std_maxj;
for sub=2:40
    odd=[odd,data(sub).v_eeg.V_EEG.odd_maxj];
    std=[std,data(sub).v_eeg.V_EEG.std_maxj];
end
odd_mean=mean(odd')';
std_mean=mean(std')';

figure
EEG_HeadTopography(data, odd_mean)     
colormap("parula"); cb=colorbar;
title(['Amplitude odd']);
ylabel(cb,'Amplitude (uV)');

figure
EEG_HeadTopography(data, std_mean)     
colormap("parula"); cb=colorbar;
title(['Amplitude std']);
ylabel(cb,'Amplitude (uV)');

figure
diff=odd_mean - std_mean;
EEG_HeadTopography(data, diff)     
colormap("parula"); cb=colorbar;
title(['Difference (odd - std) amplitude']);
ylabel(cb,'Amplitude (uV)');

% Plots related to male or female patients
m_cont=1;
f_cont=1;
m=[];
f=[];
for k=1:40
    if data(k).gender=='M'
        m(m_cont)=k;
        m_cont=m_cont+1;
    else
        f(f_cont)=k;
        f_cont=f_cont+1;
    end
end

odd_m=data(m(1)).v_eeg.V_EEG.odd_maxj;
std_m=data(m(1)).v_eeg.V_EEG.std_maxj;
for k=2:length(m)
    odd_m=[odd_m, data(m(k)).v_eeg.V_EEG.odd_maxj];
    std_m=[std_m, data(m(k)).v_eeg.V_EEG.std_maxj];
end
clear std
odd_mean_m=mean(odd_m')';
err_mean_odd_m=std(odd_m')'/sqrt(length(m));
std_mean_m=mean(std_m')';
err_mean_std_m=std(std_m')'/sqrt(length(m));

odd_f=data(f(1)).v_eeg.V_EEG.odd_maxj;
std_f=data(f(1)).v_eeg.V_EEG.std_maxj;
for k=2:length(f)
    odd_f=[odd_f, data(f(k)).v_eeg.V_EEG.odd_maxj];
    std_f=[std_f, data(f(k)).v_eeg.V_EEG.std_maxj];
end
odd_mean_f=mean(odd_f')';
err_mean_odd_f=std(odd_f')'/length(f);
std_mean_f=mean(std_f')';
err_mean_std_f=std(std_f')'/length(f);

figure
diff=odd_mean_f - std_mean_f;
EEG_HeadTopography(data, diff)     
colormap("parula"); cb=colorbar;
title(['Difference (odd - std) amplitude female subjects']);
ylabel(cb,'Amplitude (uV)');

figure
diff=odd_mean_m - std_mean_m;
EEG_HeadTopography(data, diff)     
colormap("parula"); cb=colorbar;
title(['Difference (odd - std) amplitude male subjects']);
ylabel(cb,'Amplitude (uV)');

x_odd=[odd_mean_m(n_channel),odd_mean_f(n_channel)];
error_odd=[err_mean_odd_m(n_channel),err_mean_odd_f(n_channel)];

x_std=[std_mean_m(n_channel),std_mean_f(n_channel)];
error_std=[err_mean_std_m(n_channel),err_mean_std_f(n_channel)];

xTickLabels=convertCharsToStrings(data(1).v_eeg.V_EEG.label{n_channel(1),1});
if length(n_channel)>1
    for k=2:length(n_channel)
        xTickLabels=[xTickLabels;convertCharsToStrings(data(1).v_eeg.V_EEG.label{n_channel(k),1})];
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
    errorbar(y',x_odd,error_odd,'k','linestyle','none');
    ylabel('Amplitude (uV)')
    title('Amplitude oddball events, male vs female subjects'); legend('male','female'); 
    set(gca, 'XTick', 1:length(xTickLabels), 'XTickLabel', xTickLabels);
    
subplot (2,1,2)
    b_2=bar(x_std,'grouped');
    hold on;
    [ngroups,nbars] = size(x_std);
    y = nan(nbars, ngroups);
    for i = 1:nbars
        y(i,:) = b_2(i).XEndPoints;
    end
    errorbar(y',x_std,error_std,'k','linestyle','none');
    ylabel('Amplitude (uV)')
    title('Amplitude standard events, male vs female subjects'); legend('male','female');
    set(gca, 'XTick', 1:length(xTickLabels), 'XTickLabel', xTickLabels);

p_val_odd=[];
p_val_std=[];

% Statistical test
for k=1:length(n_channel)
    dati_condizione_1 = odd_m(n_channel(k),:);  
    dati_condizione_2 = odd_f(n_channel(k),:); 
    [~, p_val_odd(k)] = ttest2(dati_condizione_1, dati_condizione_2);
end
for k=1:length(n_channel)
    dati_condizione_1 = std_m(n_channel(k),:);  
    dati_condizione_2 = std_f(n_channel(k),:); 
    [~, p_val_std(k)] = ttest2(dati_condizione_1, dati_condizione_2);

end






