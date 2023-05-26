function [model] = BPD1(m,timeseries)
%%  Geodetic Bayesian Inversion Software for Time Series (GBIS4TS) 
%   Revised by Yilin Yang, 2022
%   Institute of Earth Sciences, University of Iceland
%
%%  =========================================================================
%   BPD1 calculates the model value for a time series with one breakpoint
%   m - all the parameters for the model 
%       ('Interception'; 'Trend1'; 'TrendChange'; 'Breakpoint'; )
%   timeseries - the timeseries includes the time information
%--------------------------------------------------------------------------
time = timeseries(:,1);
time1 = time ;%- time(1);
sympref('HeavisideAtOrigin', 1);
model = m(1) - m(2)*time1(1) + m(2)*time1 + m(3)*heaviside(time-m(4)).*time ...
    - heaviside(time-m(4))*m(3).*time(find(time>=m(4),1,'first'));
%model = model - mean(model) + mean(timeseries(:,2));
end

