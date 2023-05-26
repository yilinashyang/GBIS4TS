function results = runInversion_ts(timeseries, wn_amp, invpar, model)
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
% Function that runs the MCMC Bayesian inversion for GNSS time series
% 
% Revised by Yilin Yang, Institute of Earth Sciences, University of Iceland
%
% Usage: results = runInversion(geo, gps, insar, invpar, model, modelInput, obs, nObs)
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Copyright: Marco Bagnardi, 2018
%
% Reference: 
% Bagnardi M. & Hooper A, (2018). 
% Inversion of surface deformation data for rapid estimates of source 
% parameters and uncertainties: A Bayesian approach. Geochemistry, 
% Geophysics, Geosystems, 19. https://doi.org/10.1029/2018GC007585
%
% The function may include third party software.
% =========================================================================
% Last update: 20 Sep, 2022

global outputDir  % Set global variables

%% Set starting model parameters

model.trial = model.m; % First trial using starting parameters

model.range = model.upper - model.lower; % Calculate range of m parameters

% Check that starting model is within bounds
if sum(model.m > model.upper)>0 || sum(model.m < model.lower)>0
    disp('Parameter#  Lower_Bound  Sart_Model  Upper_Bound')
    ix = find(model.m > model.upper | model.m < model.lower);
    fprintf('%d %f6 %f6 %f6', ix, model.lower(ix), model.m(ix), model.upper(ix))
    error('Starting model is out of bounds')
end

nModel = length(model.m);    % Number of model parameters to invert for
probTarget = 0.5^(1/nModel); % Target probability used in sensitivity test
probSens = zeros(nModel,1);  % Initialise vector for sensitivity test

iKeep = 0; % Initialiase kept iterations counter
iReject = 0;    % Initialise rejected iterations counter
iKeepSave = iKeep;  % Initialise saving schedule for kept iterations
iRejectSave = iReject; % Initialise saving schedule for rejected iterations

mKeep= zeros(nModel,invpar.nSave,'single');   % Initialise matrix of model parameters to keep
PKeep= zeros(1,invpar.nSave,'single');   % Initialise probability vector

POpt = -1e99; % Set initial optimal probability
%%
% Question here: Do we need simulated annealing here?
%%
T = invpar.TSchedule(1); % Set first temperature
iTemp = 0;  % Initialise temperature schedule index (for initial Simulated Annealing)
nTemp = length(invpar.TSchedule); % Number of temperature steps (for initial Simulated Annealing)
% no Simulated Annealing, so nTemp=1
%% Start core of inversion

sensitivityTest = 0; % Switch off sensitivity test at first iteration
 
while iKeep < invpar.nRuns  % While number of iterations is < than number of runs...
    if (iKeep/invpar.TRuns == round(iKeep/invpar.TRuns)) & (iTemp < nTemp) % Follow temperature schedule
        iTemp = iTemp + 1;
        T = invpar.TSchedule(iTemp);    % Assign temperature from T schedule
        % no Simulated Annealing so T=1 all the way?
        if iKeep > 0
            model.trial = model.optimal; % change parameter to optimal one (Yilin)
        end
        if T == 1                       % Set Hyperparameter when T reaches 1 (Hyperparameter is currently not is use!!!)
            setHyperParameter = 1; % always set hyperparameters?
        else
            setHyperParameter = 0;
        end
%     elseif (iKeep/invpar.TRuns == round(iKeep/invpar.TRuns)) & (iTemp == nTemp)
%         model.trial = model.optimal; % change parameter to optimal one (Yilin)
    end
    
    
    if sum(iKeep == invpar.sensitivitySchedule)>0   % Check if it's time for sensitivity test based on schedule
        sensitivityTest = 1; % Switch on sensitivity test
    end

    %% Calculate 3D displacements from model / model of Time series (Yilin) 
    nObs = length(timeseries(:,1));
    UTot = zeros(nObs);   % Initialise matrix of modeled displacements (1 x number of observation points here)
    switch invpar.model
        case 'BPD1'
            mFunc = model.trial(1:6);    % Select source model parameters from all
            U = BPD1(mFunc,timeseries);               % Calculate 3D displacements / U is the model value
        case 'BPD2'
            % implement later;
            mFunc = model.trial(1:8);    % Select source model parameters from all
            U = BPD2(mFunc,timeseries);               % Calculate 3D displacements / U is the model value
    end
       
    %% Calculate residuals
    resExp = 0; % Initialise (Gm - d)' * InvCov * (Gm - d) (Yilin: this is the squared sum of res)
    %% here I need to calculate the var-cov matrix with time-correlated noise (Yilin)
    %Cov = varMatrix(nObs, wn_amp, model.trial(end-1), model.trial(end));
    Cov = UniVarMatrix(nObs, wn_amp, model.trial(end-2), model.trial(end-1)); 
    invCov = Cov^(-1);
    rGps = timeseries(:,2) - U; % Residual
    resExp = rGps(:)' * invCov * rGps(:); % (Gm - d)' * InvCov * (Gm - d)
    
    %% Continue inversion ...
    
    if setHyperParameter == 1
        %hyperPrev = resExp/nObs; % set hyperparameter on first reaching T=1;
        hyperPrev = 1; % set hyperparameter to 1;
        model.trial(end) = log10(hyperPrev);
        setHyperParameter = 0; % Switch setHyperParameter off
    end
    
    % What is the difference between hyperPrev and hyperParam? (yilin)
    hyperParam = 1;
    
    % !! Currently hyperparameter is set to 1
    % P = -resExp/(2*hyperParam); % Probability is exp of P (logrithm)
    % here I use the eq(8) in Kouladi&Clarke(2021) to calculate P
    P = -resExp/(2) - logdet(Cov)/(2) - nObs*log(2*pi)/(2);
    
    % here I have to use eq(5) in Williams(2008) to calculate P  -- nope
    %P = -0.5*(2*nObs*log(r) + log(det(UniCov)) + resExp/(r^2) + nObs*log(2*pi));

    if iKeep>0
        PRatio = (hyperPrev/hyperParam)^(nObs/2)*exp((P-PPrev)/T);  % Probability ratio
    else
        PRatio=1; % Set to 1 for first iteration (always keep first iteration)
    end
    
    %% Perform sensitivity test if necessary and change step size
    
    if sensitivityTest > 1
        probSens(sensitivityTest-1) = PRatio; % Assign probability to current model parameter
        if sensitivityTest > nModel % Check if sensitivity test has finished
            if iKeepSave > 0
                rejectionRatio = (iReject - iRejectSave)/(iKeep - iKeepSave); % Calculate rejection rate
                probTarget = probTarget * rejectionRatio * 1/0.77; % Adjust target probability to reach 77% rejection rate
                probTarget(probTarget<1e-06) = 1e-06; % Prevent from reaching zero.
            end
            sensitivityTest = 0;    % Swtich off sensitivity test
            probSens(probSens > 1) = 1./probSens(probSens > 1);
            PDiff = probTarget - probSens;
            indexP = PDiff > 0; % Select model parameters for which to adjust model step
            model.step(indexP) = model.step(indexP).*exp(-PDiff(indexP)/probTarget*2);  % assign new model step
            indexP = PDiff < 0; % Select remaining model parameters
            model.step(indexP) = model.step(indexP).*exp(-PDiff(indexP)/(1-probTarget)*2); % assign new model step
            model.step(model.step > model.range) = model.range(model.step > model.range); % Check if step is within range
            iKeepSave = iKeep;
            iRejectSave = iReject;
        end
        
    else
        iKeep = iKeep + 1;
        if PRatio >= rand(1,1)  % If conditions are met, keep model trial
            model.m = model.trial; % Substitute m with model trial
            mKeep(:,iKeep) = model.m;   % Keep model trial
            PKeep(:,iKeep) = P;         % P of current model
            
            PPrev = P;  % Assign current P to PPrev for next trial
            hyperPrev = hyperParam; % Assign current Hyperparameter for next trial
            
            if P > POpt   % Update current optimal model if likelihood is higher
                model.optimal = model.m;
                model.funcOpt = mFunc;
                POpt = P;
            end
        else                    % Reject model trial and keep previous model
            iReject = iReject + 1;
            mKeep(:,iKeep) = mKeep(:,iKeep-1);
            PKeep(:,iKeep) = PKeep(:,iKeep-1);
        end
        
        if iKeep/invpar.nSave == round(iKeep/invpar.nSave) % display and save results at regular intervals (1000 or 10000 iterations)
            if iKeep >= 20000           % Increase time step for saving/displaying after 20000 iterations
                invpar.nSave = 10000;
            end
            
            % Print current status of inversion to screen
            disp('=========================================================')
            disp(['Model: ',invpar.model])
            disp([num2str(iKeep),' model trials. Optimal Prob = exp(',num2str(POpt),')'])
            disp(['Hyperparameter=',num2str(hyperParam)])
            disp([num2str(iReject),' models rejected: ', num2str((iReject/iKeep)*100),'% of model trials.'])
            
            % allocate space for next blocks to keep
            mKeep(:,iKeep + invpar.nSave) = 0;
            PKeep(:,iKeep + invpar.nSave) = 0;
            
            % Save results to temporary file for insepction during
            % inversion
            save([outputDir,'/temporary.mat'], 'mKeep', 'PKeep', 'model', 'timeseries', 'invpar');
            
            % Display current optimal model parameters on screen
                if invpar.model == 'BPD1'
                    fprintf('Interception: %f\n',(model.funcOpt(1,:)));
                    fprintf('Initial Trend: %f\n',(model.funcOpt(2,:)));
                    fprintf('Trend Change: %f\n',(model.funcOpt(3,:)));
                    fprintf('Breakpoint: %f\n',(model.funcOpt(4,:)));
                    fprintf('Spectral Index of noise: %f\n',(model.funcOpt(5,:)));
                    fprintf('Amplitude of Power-law noise: %f\n',(model.funcOpt(6,:)));
                elseif invpar.model == 'BPD2'
                    fprintf('Interception: %f\n',(model.funcOpt(1,:)));
                    fprintf('Initial Trend: %f\n',(model.funcOpt(2,:)));
                    fprintf('Trend Change 1: %f\n',(model.funcOpt(3,:)));
                    fprintf('Trend Change 2: %f\n',(model.funcOpt(5,:)));
                    fprintf('Breakpoint 1: %f\n',(model.funcOpt(4,:)));
                    fprintf('Breakpoint 2: %f\n',(model.funcOpt(6,:)));
                    fprintf('Spectral Index of noise: %f\n',(model.funcOpt(7,:)));
                    fprintf('Amplitude of Power-law noise: %f\n',(model.funcOpt(8,:)));
                end
        end
    end
    
    %% here I add a check for BPD2
    % for BPD1, it will continue regularly
    while 1
        if sensitivityTest > 0  % Perform sensitivity test (no models are kept during this phase!)
            randomStep = zeros(nModel,1);
            randomStep(sensitivityTest) = model.step(sensitivityTest) * sign(randn(1,1))/2; % Assign random step
            model.trial = model.m + randomStep; % New model trial
            % Check that new model trial is withing bounds
            if model.trial(sensitivityTest) > model.upper(sensitivityTest)
                model.trial(sensitivityTest) = model.trial(sensitivityTest) - model.step(sensitivityTest);
            end

            hyperParam = hyperPrev;
            sensitivityTest = sensitivityTest + 1; % Move index to that of next parameter until all parameters are done
        else
            randomStep = model.step.*(rand(nModel,1)-0.5)*2;     % Make random step
            % check whether the step for breakpoint is too small, if so, change the step to one day(yilin)
             if abs(randomStep(4)) < 0.0027
                 randomStep(4) = 0.0027 * sign(randomStep(4));
             end
            model.trial = model.m + randomStep;                 % Assign new model trial to previous + random step
            % Check that new model trial is withing bounds
            model.trial(model.trial > model.upper) = 2 * model.upper(model.trial > model.upper) - ...
                model.trial(model.trial > model.upper);

            model.trial(model.trial < model.lower) = 2 * model.lower(model.trial < model.lower) - ...
                model.trial(model.trial < model.lower);
        end
        if strcmp(invpar.model,'BPD1')
            break;
        end
        if model.trial(3) > model.trial(5)
            temp_bp = model.trial(5);
            model.trial(5) = model.trial(3);
            model.trial(3) = temp_bp;
            clear temp_bp
        end
        if abs(model.trial(5) - model.trial(3)) < (0.0027*20)
            continue;
        else
            break;
        end
    end
end

%% Clean up and prepare results
flag = find(PKeep==0);
mKeep(:, flag) = []; % Remove unused preallocated memory
PKeep(:, flag) = []; % Remove unused preallocated memory

results.mKeep = mKeep;
results.PKeep = PKeep;
results.model = model;
results.optimalmodel = model.funcOpt;

