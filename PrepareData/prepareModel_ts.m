function model = prepareModel_ts(modelCode, invpar, paraList)
%%  Geodetic Bayesian Inversion Software for Time Series (GBIS4TS) 
%   Revised by Yilin Yang, 2022
%   Institute of Earth Sciences, University of Iceland
%
%   Revised from:
%%  Geodetic Bayesian Inversion Software (GBIS)
%   Software for the Bayesian inversion of geodetic data.
%   Copyright: Marco Bagnardi, 2018
%
%%  =========================================================================
% Function to prepare model parameters
%
% Usage: model = prepareModel_ts(modelInput, invpar, timeseries)
% Input Parameters:
%       modelInput: parameters read from input file
%       invpar: inversion parameters
%
% Output Parameters:
%       model: structure containing model settings to be used for inversion
% =========================================================================
% The search range is set according to Table S2 in Yang et al.
% =========================================================================
%
% Reference: 
% Bagnardi M. & Hooper A, (2018). 
% Inversion of surface deformation data for rapid estimates of source 
% parameters and uncertainties: A Bayesian approach. Geochemistry, 
% Geophysics, Geosystems, 19. https://doi.org/10.1029/2018GC007585
%
% The function may include third party software.
% =========================================================================
%
% Revised by: Yilin Yang
%               - Institute of Earth Sciences, University of Iceland
% Last update: 7 Sep, 2022

%% Initialize model vectors
model.m = zeros(500,1);
model.step = model.m;
model.lower = model.m;
model.upper = model.m;

%% Assign model parameters from input file
% need further investigation on the custimization of the parameters
switch modelCode
    case 1 % Break point detection with one break point
        nParameters = 6;
        model.m = paraList;
        model.step =  [ 1;   paraList(2)*0.05; paraList(2)*0.05; 0.0027; 0.05; 1];
        model.lower = [-5;   paraList(2)*(-1); paraList(3)*(-1); paraList(4)-1; -1.5;               0];
        model.upper = [ 5;   paraList(2)*( 2); paraList(3)*( 2); paraList(4)+1;    0; paraList(6)*1.5];
        % for negative velocity, change the boundary
        if model.lower(2) > model.upper(2)
            temp = model.lower(2);
            model.lower(2) = model.upper(2);
            model.upper(2) = temp;
        end
        % for negative velocity, change the boundary
        if model.lower(3) > model.upper(3)
            temp = model.lower(3);
            model.lower(3) = model.upper(3);
            model.upper(3) = temp;
        end
        model.gaussPrior = false(nParameters,1);
        model.parName = {'Intercept';'Trend1'; 'TrendChange'; ...
                         'Breakpoint'; 'SpectralIndex'; 'Amplitude'};

    case 2 % Break point detection with two break points
        nParameters = 8;
        model.m = paraList;
        model.step = [  1;paraList(2)*0.05; paraList(2)*0.05; 0.0027; paraList(2)*0.05; 0.0027; 0.05; 1];
        model.lower = [-5;paraList(2)*(-1); paraList(3)*(-1); paraList(4)-0.5; ...
                       paraList(5)*(-1); paraList(6)-0.5; -1.5; 0];
        model.upper = [ 5;paraList(2)*( 2); paraList(3)*( 2); paraList(4)+0.5; ...
                       paraList(5)*( 2); paraList(6)+0.5; 0; paraList(8)*1.5];
        % for negative velocity, change the boundary
        if model.lower(2) > model.upper(2)
            temp = model.lower(2);
            model.lower(2) = model.upper(2);
            model.upper(2) = temp;
        end
        % for negative velocity, change the boundary
        if model.lower(3) > model.upper(3)
            temp = model.lower(3);
            model.lower(3) = model.upper(3);
            model.upper(3) = temp;
        end
        % for negative velocity, change the boundary
        if model.lower(4) > model.upper(4)
            temp = model.lower(4);
            model.lower(4) = model.upper(4);
            model.upper(4) = temp;
        end
        % for negative velocity, change the boundary
        if model.lower(5) > model.upper(5)
            temp = model.lower(5);
            model.lower(5) = model.upper(5);
            model.upper(5) = temp;
        end
        if model.lower(6) > model.upper(6)
            temp = model.lower(6);
            model.lower(6) = model.upper(6);
            model.upper(6) = temp;
        end
        model.gaussPrior = false(nParameters,1);
        model.parName = {'Interception'; 'Trend1'; 'TrendChange1'; ...
                         'Breakpoint1'; 'TrendChange2'; 'Breakpoint2'; ...
                         'SpectralIndex'; 'Amplitude'};

        % --------------------Comments by Yilin----------------------------
        % Removed all other models here, in case of errors or confilcts
        % Other time series models can be added here
        % e.g. for time-variable periodc signals
        % -----------------------------------------------------------------
    otherwise
        error('Invalid model')
end

% there´s a hyperparameter in the end of the vector
model.m(nParameters+1) = 0;
model.step(nParameters+1) = 1e-3;
model.lower(nParameters+1) = -5e-1;
model.upper(nParameters+1) = 5e-1;

% Clear unused values
model.m = model.m(1:nParameters+1);
model.step = model.step(1:nParameters+1);
model.lower = model.lower(1:nParameters+1);
model.upper = model.upper(1:nParameters+1);

