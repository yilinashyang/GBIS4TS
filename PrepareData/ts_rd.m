function [ Timeseries,un ] = ts_rd( filename )
%%  Geodetic Bayesian Inversion Software for Time Series (GBIS4TS) 
%   by Yilin Yang, 2022
%   Institute of Earth Sciences, University of Iceland
%
%%  =======================================================================
% This Function is used to input a txt format file and output a timeserise
%
% Updated on 2 March 2023
%%
[YD,N,un]=textread(filename,'%f %f %f','headerlines',1); %#ok<*DTXTRD>
Timeseries=[YD,N];
end

