function [cMA,cAR]=rico_mod(z,B,fn,fc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters
%	z  (0.01) minimum attenuation at frequency fc
%	B  (2-8) bandwidth corresponding to attenuation 0.707
% 	fn (50) noise frequency
% 	fc  sampling frequency
% 
% Output parameters
%	cMA  filter coefficients (MA part)
%	cAR  filter coefficients (AR part)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=1/fc;
b=pi*B*T;
a=b*z;
c1=-2*(1-a)*cos(2*pi*fn*T);
c2=(1-a)^2;
c3=2*(1-b)*cos(2*pi*fn*T);
c4=-(1-b)^2;
cMA=[1 c1 c2];
cAR=[1 -c3 -c4];

end