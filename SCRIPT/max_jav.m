function    data=max_jav(data,sub,fc,epoch_camp,n_before)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for Maximum Amplitude Calculation with Jitter Correction
%
% Input:        data: Struct containing EEG data
%               sub: Subject index
%               fc: Sampling frequency
%               epoch_camp: Number of samples in one epoch
%               n_before: Number of samples before event for baseline correction
%
% Output:       data: Updated struct with maximum amplitude and error for oddball and standard events
%
% Authors:      Roberto Pilotto
%               Salvatore Rapisarda
%               Chiara Razzuoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jitter is defined as the difference between all epoch latencies of the channel and the latency of the first epoch
delay_odd=round(data(sub).v_eeg.V_EEG.epoch_odd_latency*fc);
for k=1:30
    delay_odd(k,:)=delay_odd(k,:)-delay_odd(k,1);
end

% ODDBALL
% Epochs division into three subgroups, each containing 1/3 of the analysed events
l_odd=round( length( delay_odd(k,:) )/3 );
for k=1:30
    [epoch_odd1(k,:),err_odd1(k,:)]=av_jitter(data(sub).v_eeg.V_EEG.epoch_odd(k,1:l_odd*epoch_camp),l_odd,epoch_camp,delay_odd(k,1:l_odd));  
end
% Baseline correction
base_corr=mean(epoch_odd1(:,1:n_before)')';
% Amplitude research
[max_odd1,ind]=max(epoch_odd1');
max_odd1=max_odd1'-base_corr;
for k=1:30
    max_err_odd1(k,1)=err_odd1(k,ind(k));
end

for k=1:30
    [epoch_odd2(k,:),err_odd2(k,:)]=av_jitter(data(sub).v_eeg.V_EEG.epoch_odd(k,l_odd*epoch_camp+1:2*l_odd*epoch_camp),l_odd,epoch_camp,delay_odd(k,l_odd+1:2*l_odd));  
end
% Baseline correction
base_corr=mean(epoch_odd2(:,1:n_before)')';
% Amplitude research
[max_odd2,ind]=max(epoch_odd2');
max_odd2=max_odd2'-base_corr;
for k=1:30
    max_err_odd2(k,1)=err_odd2(k,ind(k));
end

l_odd_provv=length(delay_odd(1,(2*l_odd+1):end));
if l_odd_provv~=1 
    for k=1:30
        [epoch_odd3(k,:),err_odd3(k,:)]=av_jitter( data(sub).v_eeg.V_EEG.epoch_odd( k,(2*l_odd*epoch_camp+1):end ),l_odd_provv,epoch_camp,delay_odd( k,(2*l_odd+1):end ));  
    end
else
    epoch_odd3=data(sub).v_eeg.V_EEG.epoch_odd( :,(2*l_odd*epoch_camp+1):end );
end
% Baseline correction
base_corr=mean(epoch_odd3(:,1:n_before)')';
% Amplitude research
[max_odd3,ind]=max(epoch_odd3');
max_odd3=max_odd3'-base_corr;
if l_odd_provv~=1 
    for k=1:30
        max_err_odd3(k,1)=err_odd3(k,ind(k));
    end
else
    max_err_odd3=zeros(30,1);
end



data(sub).v_eeg.V_EEG.odd_max3j=[max_odd1,max_odd2,max_odd3];
data(sub).v_eeg.V_EEG.odd_max3j_err=[max_err_odd1,max_err_odd2,max_err_odd3];
% Mean of the 3 subgroups
data(sub).v_eeg.V_EEG.odd_maxj=mean(data(sub).v_eeg.V_EEG.odd_max3j')';
% Standard error
data(sub).v_eeg.V_EEG.odd_maxj_err=std(data(sub).v_eeg.V_EEG.odd_max3j')'/sqrt(3);


% STANDARD
% Epochs division into three subgroups, each containing 1/3 of the analysed events
delay_std=round(data(sub).v_eeg.V_EEG.epoch_std_latency*fc);
for k=1:30
    delay_std(k,:)=delay_std(k,:)-delay_std(k,1);
end
% delay_std=int(delay_std);
l_std=round( length( delay_std(k,:) )/3 );
for k=1:30
    [epoch_std1(k,:),err_std1(k,:)]=av_jitter(data(sub).v_eeg.V_EEG.epoch_std(k,1:l_std*epoch_camp),l_std,epoch_camp,delay_std(k,1:l_std));  
end
% Baseline correction
base_corr=mean(epoch_std1(:,1:n_before)')';
% Amplitude research
[max_std1,ind]=max(epoch_std1');
max_std1=max_std1'-base_corr;
for k=1:30
    max_err_std1(k,1)=err_std1(k,ind(k));
end

for k=1:30
    [epoch_std2(k,:),err_std2(k,:)]=av_jitter(data(sub).v_eeg.V_EEG.epoch_std(k,l_std*epoch_camp+1:2*l_std*epoch_camp),l_std,epoch_camp,delay_std(k,l_std+1:2*l_std));  
end
% Baseline correction
base_corr=mean(epoch_std2(:,1:n_before)')';
% Amplitude research
[max_std2,ind]=max(epoch_std2');
max_std2=max_std2'-base_corr;
for k=1:30
    max_err_std2(k,1)=err_std2(k,ind(k));
end

l_std_provv=length(delay_std(1,(2*l_std+1):end));
if l_std_provv~=1 
    for k=1:30
        [epoch_std3(k,:),err_std3(k,:)]=av_jitter( data(sub).v_eeg.V_EEG.epoch_std( k,(2*l_std*epoch_camp+1):end ),l_std_provv,epoch_camp,delay_std( k,(2*l_std+1):end ));  
    end
else
    epoch_std3=data(sub).v_eeg.V_EEG.epoch_std( :,(2*l_std*epoch_camp+1):end );
end
% Baseline correction
base_corr=mean(epoch_std3(:,1:n_before)')';
% Amplitude research
[max_std3,ind]=max(epoch_std3');
max_std3=max_std3'-base_corr;
if l_std_provv~=1 
    for k=1:30
        max_err_std3(k,1)=err_std3(k,ind(k));
    end
else
    max_err_std3=zeros(30,1);
end



data(sub).v_eeg.V_EEG.std_max3j=[max_std1,max_std2,max_std3];
data(sub).v_eeg.V_EEG.std_max3j_err=[max_err_std1,max_err_std2,max_err_std3];
% Mean of the 3 subgroups
data(sub).v_eeg.V_EEG.std_maxj=mean(data(sub).v_eeg.V_EEG.std_max3j')';
% Standard error
data(sub).v_eeg.V_EEG.std_maxj_err=std(data(sub).v_eeg.V_EEG.std_max3j')'/sqrt(3);

end
