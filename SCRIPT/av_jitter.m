function [y_aj,err_aj]=av_jitter(x,n_ep,camp_ep,ritardi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to apply the averaging technique with jitter correction
%
% Input:    x:          Signal to be filtered
%           n_ep:       Number of epochs to average
%           camp_ep:    Number of samples in one epoch
%           ritardi:    Vector containing delays from the second
%                       to the last epoch relative to the first
%
% Output:   y_aj:       Filtered signal with jitter correction
%           err_aj:     Standard error of the averaged epoch
%
% Authors:  Roberto Pilotto
%           Salvatore Rapisarda
%           Chiara Razzuoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adjust the input data dimensions to make the function as general as possible
[righe,colonne]=size(x);
if (righe==n_ep)
    y=x;
elseif (colonne==n_ep)
    y=x';
elseif (colonne==1)
    y=reshape(x,[camp_ep,n_ep])';
elseif (righe==1)
    y=reshape(x,[camp_ep,n_ep])';
end

y_j=y;
for k=1:n_ep-1
    % jitter correction
    y_j(k+1,:)=circshift(y(k+1,:),-ritardi(k));
end
err_aj=std(y_j)/sqrt(n_ep);
y_aj=mean(y_j);
end