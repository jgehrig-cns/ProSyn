function func_TPSim_synPros(ID, tlock)

params             = [];
params.cndCmb      = {'TT', 'TT';
 'IT', 'IT';
 'TN', 'TN'}; 

% {'TT', 'TN';
%  'TT', 'IT';
%  'TN', 'IT'}; % N x 2 cell array with 'TT', 'TN', 'IT'. 
  % If 'statistics' = 'yes', specify two different conditions per row

tlock_options = {'stimOnset',  'stimOffset', 'speechOnset',
                 [-2:0.01:6], [-4:0.01:5],  [-5.6:0.01:2.4]};
   
params.timelock    = tlock_options(1,tlock);     % determines which TFR data are loaded. Can be 'offset'
params.time        = tlock_options{2,tlock}; % -2:0.01:6;    % begin : step : end



% params.timelock    = 'speech_onset'; %'offset'     % determines which TFR data are loaded. Can be 'offset'
params.method      = 'TPSim';        % 'TPSim', 'corrDist', 'euclDist'
params.subjects    = ID;              % subject ID or 'all'
params.freq        = 1:48;           % freq1:freq2
% params.time        = -5.6:0.01:2.4;  %    -4:0.01:5;    % begin : step : end %-4:0.01:5
params.cmbExhaust  = 'no';           % 'yes' or 'no', if 'yes', the calculations are performed
% on (cond1, cond2) x (cond1, cond2) and not only
% on cond1 x cond2. Only applicable if cond1 ~= cond2.
params.winSize     = 0.5;            % window size at 4Hz if taper = 'yes' and at all freqs if taper = 'no'
params.stepWidth   = 0.1;          % time step width (movement of window per step) at 4Hz if taper = 'yes' and at all freqs if taper = 'no'
params.winTaper    = 'no';         % if 'linear', tapering of window size according to $ winSize*freq(1)/freq $ is applied
params.stepTaper   = 'no';         % if 'linear', tapering of time step according to $ stepWidth*freq(1)/freq $ is applied
params.winLags     = 0;            % array of window lags. Should be symmetric
% around 0, because only one condition is lagged.
% [-0.01, -0.005, 0, 0.005, 0.01];

% METHOD-SPECIFIC PARAMETERS
%---------------------------

% TPSim:
params.absolute    = 'no';         % 'yes' if mean over trial combinations
% should be taken over absolute values,
% 'no' if signs should be included


% Set the subjects
if ischar(params.subjects) && strcmp(params.subjects, 'all')
    subjects = 1:9;
elseif isnumeric(params.subjects)
    subjects = params.subjects;
else
    error('No such subject.');
end

%------------------------------------------------
% Set the MATLAB path
%------------------------------------------------

restoredefaultpath
% load and start the fieldtrip toolbox
addpath(strcat('/data', ...
    '/projects', ...
    '/ecog', ...
    '/fieldtrip', ...
    '/fieldtrip-20210328'));
ft_defaults;
% add path with design file
addpath(strcat('/data', ...
    '/projects', ...
    '/ecog', ...
    '/skript', ...
    '/new_stim'));
% add path to external functions
addpath(strcat('/data', ...
    '/projects', ...
    '/ecog', ...
    '/Mastermodul', ...
    '/functions', ...
    '/external'));
% add path to TPSim functions
addpath(strcat('/data', ...
    '/projects', ...
    '/ecog', ...
    '/Mastermodul', ...
    '/functions', ...
    '/TPSim'));
% add path to TFR files
addpath(strcat('/data', ...
    '/projects', ...
    '/ecog', ...
    '/Mastermodul', ...
    '/files', ...
    '/TFR', ...
    '/2024_03_05_TFRwaveletFourier'));
% add path to TPSim files
addpath(strcat('/data', ...
    '/projects', ...
    '/ecog', ...
    '/Mastermodul', ...
    '/files', ...
    '/TPSim', ...
    '/2021_08_11_TPSimFromWavelet180'));

%------------------------------------------------
% Load necessary files
%------------------------------------------------

% Load design file
load(strcat('/data', ...
    '/projects', ...
    '/ecog', ...
    '/skript', ...
    '/new_stim', ...
    '/design2021_2.mat'));

%------------------------------------------------
% Loop over subjects
%------------------------------------------------

for subject = subjects
    % Load time-frequency data
        loadname = (strcat('/data', ...
            '/projects', ...
            '/ecog', ...
            '/Mastermodul', ...
            '/files', ...
            '/TFR', ...
            '/2024_03_05_TFRwaveletFourierSynPros/TFRwaveletFourier_', ...
            Design{subject, 3}, '_', ...
            params.timelock,'.mat'));
  load(loadname{:})
    % Store TFR data
    TFR = TFRwave; clear TFRwave
    
    %------------------------------------------------
    % Loop over condition combinations
    %------------------------------------------------
    
    % Find number of combinations
    [nCombs, ~] = size(params.cndCmb);
    % prepare temporary structure
    tmp = [];
    for cmb = 1:nCombs
        % derive condition names of current combination
        if nCombs == 1
            tmp.cnd1 = params.cndCmb{1};
            tmp.cnd2 = params.cndCmb{2};
        else
            tmp.cnd1 = params.cndCmb{cmb, 1};
            tmp.cnd2 = params.cndCmb{cmb, 2};
        end
        % the time of interest
        tmp.time  = params.time;
        % the channel labels
        tmp.label = TFR.data.label;
        % replace condition names with indices 1='TT', 2='TN', 3='IT'
        if strcmp(tmp.cnd1, 'TT'); tmp.c1 = 1; end
        if strcmp(tmp.cnd1, 'TN'); tmp.c1 = 2; end
        if strcmp(tmp.cnd1, 'IT'); tmp.c1 = 3; end
        if strcmp(tmp.cnd2, 'TT'); tmp.c2 = 1; end
        if strcmp(tmp.cnd2, 'TN'); tmp.c2 = 2; end
        if strcmp(tmp.cnd2, 'IT'); tmp.c2 = 3; end
        % find trials of both conditions
        tmp.c1 = find(condTrials == tmp.c1);
        tmp.c2 = find(condTrials == tmp.c2);
        
        %------------------------------------------------
        % Loop over frequencies
        %------------------------------------------------
        
        for freq = params.freq
            disp(freq)
            % Set frequency of interest
            tmp.freq = freq;
            % the window size tapering
            if strcmp(params.winTaper, 'linear')
                tmp.winSize = params.winSize*(params.freq(1)/tmp.freq);
            elseif strcmp(params.winTaper, 'no')
                tmp.winSize = params.winSize;
            else
                error('Window taper has to be either "linear" or "no".');
            end
            % the time step tapering
            if strcmp(params.stepTaper, 'linear')
                tmp.stepWidth = params.stepWidth*(params.freq(1)/tmp.freq);
            elseif strcmp(params.winTaper, 'no')
                tmp.stepWidth = params.stepWidth;
            else
                error('Time step taper has to be either "linear" or "no".');
            end
            % the time points at which the windows should start
            tmp.tSteps = tmp.time(1) : tmp.stepWidth : (tmp.time(end) - tmp.winSize);
            % initialise container for results and lags
            if tmp.freq == params.freq(1)
                tmp.datAll = zeros(length(tmp.label), ...
                    length(tmp.tSteps), ...
                    params.freq(end));
                tmp.lagAll = tmp.datAll;
            end
            
            %------------------------------------------------
            % loop over channels
            %------------------------------------------------
            
            for chan = 1:length(tmp.label)
                % store the channel label
                tmp.chan      = chan;
                % counter for the time step
                tmp.nSteps    = 0;
                % find index of frequency within TFR
                [~, tmp.fInd] = min(abs(TFR.data.freq - tmp.freq));
                % TFR dimord = 'rpttap_channel_freq_time'
                tmp.TFR       = squeeze(TFR.data.fourierspctrm(:, ...
                    tmp.chan, ...
                    tmp.fInd, ...
                    :));
                
                %------------------------------------------------
                % loop over time steps
                %------------------------------------------------
                
                for t = tmp.tSteps
                    % Save the current time
                    tmp.t      = t;
                    % increase the time step counter
                    tmp.nSteps = tmp.nSteps + 1;
                    % initialise container for result and lag
                    tmp.result = 0;
                    tmp.maxLag = 0;
                    
                    %------------------------------------------------
                    % loop over time lags
                    %------------------------------------------------
                    
                    for lag = params.winLags
                        % store the lag
                        tmp.lag        = lag;
                        % Set temporal boundary indices for first condition's time window
                        [~, tmp.t1]    = min(abs(TFR.data.time -  tmp.t));
                        [~, tmp.t2]    = min(abs(TFR.data.time - (tmp.t + tmp.winSize)));
                        % Set temporal boundary indices for second condition's(lagged) time window
                        [~, tmp.t1Lag] = min(abs(TFR.data.time - (tmp.t + tmp.lag)));
                        [~, tmp.t2Lag] = min(abs(TFR.data.time - (tmp.t + tmp.lag + tmp.winSize)));
                        % adjust different lengths of the time windows. Take the
                        % minimum of timesteps as common length
                        tmp.minSteps   = min([length(tmp.t1:tmp.t2), length(tmp.t1Lag:tmp.t2Lag)]);
                        % choose the same number of steps from both sets
                        tmp.stepsFix   = tmp.t1:tmp.t2;
                        tmp.stepsFix   = tmp.stepsFix(1:tmp.minSteps);
                        tmp.stepsLag   = tmp.t1Lag:tmp.t2Lag;
                        tmp.stepsLag   = tmp.stepsLag(1:tmp.minSteps);
                        % Store the TFR data of the condition's trials in two matrices
                        if strcmp(params.cmbExhaust, 'no')
                            tmp.tr1    = tmp.c1;
                            tmp.tr2    = tmp.c2;
                        elseif strcmp(params.cmbExhaust, 'yes') && ~strcmp(tmp.cnd1, tmp.cnd2)
                            tmp.tr1    = [tmp.c1; tmp.c2];
                            tmp.tr2    = [tmp.c1; tmp.c2];
                        else
                            error('cmbExhaust not specified correctly.');
                        end
                        % load the TFR data of the chosen trials and the time steps
                        A              = tmp.TFR(tmp.tr1, tmp.stepsFix);
                        B              = tmp.TFR(tmp.tr2, tmp.stepsLag);
                        % container for the result
                        res            = zeros(length(A(:, 1)), ...
                            length(B(:, 1)));
                        
                        %------------------------------------------------
                        % loop over trial combinations
                        %------------------------------------------------
                        
                        if strcmp(params.method, 'TPSim')
                            % convolve each tuple of trials at each timestep in
                            % the time window and take the mean thereof
                            for i = 1:length(A(:, 1))
                                for j = 1:length(B(:, 1))
                                    % Define the correlation matrix according to
                                    %
                                    %                   tr1(n)  x     tr2(m)*
                                    %   corr(n,m) = --------------------------
                                    %               abs(tr1(n)) x abs(tr2(m)*)
                                    %
                                    % where * denotes the complex conjugate.
                                    res(i, j) = nanmean(...
                                        conj(A(i, :))./abs(conj(A(i, :))).*...
                                        B(j, :)./abs(B(j, :)));
                                end
                            end
                            % Only take the real part
                            res = real(res);
                            % If wanted, take the absolute value before
                            % averaging
                            if strcmp(params.absolute, 'yes')
                                res = abs(res);
                            end
                            % fisher-transform the data
                            % "1's" on the diagonal are removed on the way
                            res = fisherz(res(res < 0.99));
                            % take the mean
                            res = squeeze(nanmean(res));
                        end
                        if strcmp(params.method, 'corrDist')
                            % convolve each tuple of trials at each timestep in
                            % the time window and take the mean thereof
                            for i = 1:length(A(:, 1))
                                for j = 1:length(B(:, 1))
                                    % Define the correlation matrix according to
                                    %
                                    %                   tr1(n)  x     tr2(m)*
                                    %   corr(n,m) = --------------------------
                                    %               abs(tr1(n)) x abs(tr2(m)*)
                                    %
                                    % where * denotes the complex conjugate.
                                    res(i, j) = nanmean(...
                                        conj(A(i, :))./abs(conj(A(i, :))).*...
                                        B(j, :)./abs(B(j, :)));
                                end
                            end
                            % Only take the real part
                            res = real(res);
                            % fisher-transform the data
                            % "1's" on the diagonal are removed on the way
                            res = fisherz(res(res < 0.99));
                            % take the mean of the absolute value of TPSim
                            % subtracted from 1
                            res = squeeze(nanmean(1 - abs(res)));
                            disp(res)
                        end
                        % check, if result is larger with the actual lag than with
                        % previous ones
                        if abs(res) > abs(tmp.result)
                            tmp.result = res;
                            tmp.maxLag = tmp.lag;
                        end
                    end
                    % Save the actual result and the corresponding lag
                    tmp.datAll(tmp.chan, tmp.nSteps, tmp.freq) = tmp.result;
                    tmp.lagAll(tmp.chan, tmp.nSteps, tmp.freq) = tmp.maxLag;
                end
            end
        end
        clear A B tmp.stepsFix tmp.stepsLag tmp.minSteps
        clear tmp.t1 tmp.t2 tmp.t1Lag tmp.t2Lag tmp.tr1 tmp.tr2
        clear tmp.lag fInd tmp.chan tmp.nSteps tmp.TFR
        clear tmp.result tmp.maxLag tmp.t tmp.c1 tmp.c2
        
        %------------------------------------------------
        % Save the results
        %------------------------------------------------
        
        if strcmp(params.method, 'TPSim')
            TPSim        = [];
            TPSim.data   = tmp.datAll;
            TPSim.info   = tmp;
            TPSim.params = params;
            TPSim.lag    = tmp.lagAll;
            % set the method label
            methLabel    = 'TPSim';
        end
        if strcmp(params.method, 'corrDist')
            TPSim        = [];
            TPSim.data   = tmp.datAll;
            TPSim.info   = tmp;
            TPSim.params = params;
            TPSim.lag    = tmp.lagAll;
            % set the method label
            methLabel    = 'corrDist';
        end
        % set the taper label
        if strcmp(params.winTaper, 'linear')
            taperLabel = '_linTaper';
        else
            taperLabel = '';
        end
        % set the lag label
        if params.winLags == 0
            lagLabel   = '';
        else
            lagLabel   = '_lag';
        end
        % set the condition label
        if strcmp(tmp.cnd1, tmp.cnd2)
            condLabel  = tmp.cnd1;
        else
            condLabel  = strcat(tmp.cnd1, '_', tmp.cnd2);
        end
        % set the timelock label
        tLockLabel = params.timelock;
        % save the results
        savename = strcat('/data', ...
            '/projects', ...
            '/ecog', ...
            '/Mastermodul', ...
            '/files', ...
            '/TPSim', ...
            '/2024_02_28_TPSimFromWavelet180_SynPros', ...
            '/', ...
            methLabel, '_', ...
            condLabel, '_', ...
            Design{subject, 3}, '_', ...
            tLockLabel, ...
            '.mat');
      save(savename{:},'TPSim');
    end
    clear tmp
end
clear subject

