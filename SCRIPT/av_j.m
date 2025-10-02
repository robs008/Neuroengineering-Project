function [y_a]=av_j(x,n_ep,camp_ep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to apply the averaging technique 
%
% Input:    x:          Signal to be filtered
%           n_ep:       Number of epochs to average
%           camp_ep:    Number of samples in one epoch
%
% Output:   y_a:        Averaged signal
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

y_a=mean(y);

end