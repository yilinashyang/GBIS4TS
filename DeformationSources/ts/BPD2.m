function [model] = BPD2(m,timeseries)
%%  Geodetic Bayesian Inversion Software for Time Series (GBIS4TS) 
%   Revised by Yilin Yang, 2022
%   Institute of Earth Sciences, University of Iceland
%
%%  =========================================================================
%   BPD2 calculates the model value for a time series with two breakpoint
%   m - all the parameters for the model 
%       ('Interception' 'Trend1'; 'TrendChange1'; 'Breakpoint1'; 'TrendChange2'; 'Breakpoint2';)
%   timeseries - the timeseries includes the time information
%--------------------------------------------------------------------------
time = timeseries(:,1);
time1 = time;% - time(1);
sympref('HeavisideAtOrigin', 1);
model = m(1) - m(2)*time1(1) + m(2).*time1 + m(3)*heaviside(time1-m(4)).*time1 ...
    - heaviside(time1-m(4))*m(3).*time1(find(time1>=m(4),1,'first')) ...
    + m(5)*heaviside(time1-m(6)).*time1 - heaviside(time1-m(6))*m(5).*time1(find(time1>=m(6),1,'first'));
%model = model - mean(model) + mean(timeseries(:,2));
end
