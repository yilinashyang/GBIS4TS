function [var_cor] = UniVarMatrix(m, wn_amp, kappa, pln_amp)
%%  Geodetic Bayesian Inversion Software for Time Series (GBIS4TS) 
%   by Yilin Yang, 2022
%   Institute of Earth Sciences, University of Iceland
%
%%  =======================================================================
%   VARMATRIX generates the variance-covariance matrix for time series 
%   Revised by Yilin Yang, at University of Iceland.
%
%   Reference: 
%       Williams S.D.P. 2003. The effect of coloured noise on the uncertainties
%       of rates estimated from geodetic time series. Journal of Geodesy. 
%   Last updated: 
%       Sep 21, 2022
%%  =======================================================================
% this is to create the time-correlation
b=zeros(m,1);
kk=kappa/2;
b(1) = 1;
for i=2:m
    b(i)=((i-1-1-kk)/(i-1))*b(i-1);
end
for i=1:m-1
    T(i+1:m,i+1)=b(1:m-i);
end
T(1:m,1)=b;
T1=(1/365)^(-1*kappa/4)*T; % T1 is the trasformation matrix
% this section is abandon as it causing infinite value of det(Cov)
% the cov matrix for the powerlaw
Cp = T1*T1';
% add all together (eq 4 in Williams 2003 JG paper)
var_cor = wn_amp*wn_amp*eye(m) + pln_amp*pln_amp*Cp;

% %% here we apply the solution from Williams S.D.P 2008. GPS Solut. (CATS paper)
% %%more detail see eq (4) - (8)
% % calculatet the cov matrixes
% J = T1*T1'; % for the power law
% I = eye(m); % for the white noise
% 
% % calculate the angle
% phi = atan2(pln_amp, wn_amp);
% 
% % calculate the scalar r
% sigma = wn_amp/cos(phi);
% 
% % calculate the 1-D covariance matrix
% UniCovMatrix = cos(phi)^2*I + sin(phi)^2*J;

end

