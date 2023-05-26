function [] = GBISrun_ts(timeseriesList, startParameter, WNlist, modelCode, nRuns)
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
%   Usage: GBISrun(timeseriesList, startParameter, modelCode, nRuns, Direction)
%
%   timeseriesList:  a file with a list of time series files
%                    format: time(yr) Coordinate Uncertainty (in .txt)
%
%   startParameter:  a file with a list of start parameters for each time series
%                    format: site_name intercep v dv epoch kappa amp
%                      unit:     // mm mm mm decimal_year // mm/(yr)^x
%                                                       x depends on kappa
%                    example: KASC 2020 21 -5 2020.0000 -1 8
%
%   WNlist: a file with the white noise amplitude for each station
%           !!Please keep the same order as timeseries list!!
%
%   modelCode:      select models to use;
%                   '1' for breakpoint detection with 1 break point
%                   '2' for breakpoint detection with 2 break point
%
%   nRuns:          number of iterations (samples) to be performed (e.g., 1000000)
%   
%   =========================================================================
%
%   Reference:
%   Bagnardi M. & Hooper A, (2018).
%   Inversion of surface deformation data for rapid estimates of source
%   parameters and uncertainties: A Bayesian approach. Geochemistry,
%   Geophysics, Geosystems, 19. https://doi.org/10.1029/2018GC007585
%
%   Update on:  20 Sep, 2022
%   =========================================================================
%% Check number of input arguments and return error if not sufficient

if nargin == 0
    help GBISrun;
    return;
end

if nargin < 4
    disp('#############################################')
    disp('####### Not enough input arguments. #########')
    disp('#############################################')
    disp(' ')
    disp('Type "help GBISrun" for more information');
    disp(' ')
    return;
end

%% Start timer and initialise global variables
clc % Clean screen
tic % Start timer

clear  outputDir  % Clear global variables
global outputDir  % Set global variables

% I don't print header here (Yilin) - maybe add later.

%% Read time series path (Yilin)
% I suggest to have all the file named as 'XXXX.txt' XXXX is the site name.
tsList = textread(timeseriesList,'%s'); % cell matrix with all the files
%% Read the fixed WN amplitude (Yilin) (unit: mm)
% Please make sure the order of WN is the same with tsList
wn = textread(WNlist,'%f'); % cell matrix with all the files

%% Read start parameters (Yilin)
% I don't explain every variable here as they are more or less obvious
% e.g. intList is the list of interception
% the number of parameters depends on the model (number of bp)
% Please KEEP THE SAME ORDER with tsList
if modelCode == 1
    [siteList,intList,vList,dvList,epoList,kList,ampList]...
        =textread(startParameter,'%s\t%f\t%f\t%f\t%f\t%f\t%f');
elseif modelCode == 2
    [siteList,intList, vList,dv1List,epo1List,dV2List,epo2List,kList,ampList]...
        =textread(startParameter,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f');
end

%% Create output directories 
outputDir = pwd; % I set the output folder as current folder
saveName = ['Output' num2str(modelCode)];
outputDir = [outputDir,'/',saveName];
disp(['Output directory: ', outputDir])
mkdir(outputDir)

%% Loop for all the stations  (Yilin)
for i = 1:length(tsList)
    %% Initialise variables (Yilin)
    timeseries = []; % Time series
    paraList = [];
    nObs = 0;   % Initialise number of observations variable
    % Obtain the site name, usually the first four digits of the file name
    currentSite = char(tsList(i)); currentSite = currentSite(end-7:end-4);
    if ~strcmp(currentSite, char(siteList(i)))
        fprintf('The starPara does not match with the ts list!\n');return;
    end  
    %% Read corresponding time series (Yilin)
    % I don't use the uncertainties in the time series because sometimes it
    % makes less sense.
    [timeseries, ~] = ts_rd(char(tsList(i)));
    %timeseries(:,2) = timeseries(:,2) - timeseries(1,2); % force the intercept to be zero.
    nObs = length(timeseries(:,1)); % calculate the epoch number

    %% Setup inversion parameters
    fprintf('Preparing estimation for %s...\n', currentSite);
    
    %??
    invpar.TRuns = 1000; % Number of runs at each temperature (Simulated Annealing only)

    invpar.nSave = 1000;    % Save output to file every 1000 iterations (every 10,000 after 20,000 iterations)
    invpar.sensitivitySchedule = [1:100:10000,11000:1000:30000,40000:10000:nRuns]; % sensitivity schedule (when to change step sizes)

    % I skip simulated annealing for now. add later if needed. 
    %invpar.TSchedule = 1; % No temperature schedule if Simulated Annealing is not performed
    % here I add the simulated annealing back for test
    invpar.TSchedule = 10.^(3:-0.2:0); % Cooling schedule for Simulated Annealing
    %
    invpar.nRuns = nRuns;
    if modelCode == 1
        invpar.model='BPD1';
        paraList = [intList(i);vList(i);dvList(i);epoList(i);kList(i);ampList(i)];
    elseif modelCode == 2
        invpar.model='BPD2';
        paraList = [intList(i);vList(i);dv1List(i);epo1List(i);...
            dV2List(i);epo2List(i);kList(i);ampList(i);];
    end
    model = prepareModel_ts(modelCode, invpar, paraList);
    %% Run inversion

    invResults = runInversion_ts(timeseries, wn(i), invpar, model);
    %% Create *.mat file with final results

    cd(outputDir)
    save([char(siteList(i,:)) '.mat'], ...
        'model', 'timeseries', 'invpar', 'invResults')
    delete temporary.mat
    cd ..

end
end