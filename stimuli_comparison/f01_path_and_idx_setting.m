function [paths,idx, ps] = f01_path_and_idx_setting(root_path)
%% paths ... 
paths.stimuli = fullfile(root_path,'f01_stimuli'); 

results_root_path = 'f02_stimuli_comparison_results';

paths.audio_spectrogram = fullfile(root_path,results_root_path,'r01_stimuli_spectrogram');
making_folder(paths.audio_spectrogram);

paths.envelopes = fullfile(root_path,results_root_path,'r02_envelopes');
making_folder(paths.envelopes);

paths.env_freq_bayesian = fullfile(paths.envelopes,'bayesian_env_fft');
making_folder(paths.env_freq_bayesian);


%% idx ... 
idx.plot_spectrogram = 0;%# % plot spectrogram ... 
idx.envelope_info_extraction = 1;%# % ... 
idx.envelope_freq_statistics = 0;%# % ...
idx.RSM = 0;%# % ... 
idx.time_domain_statistics = 0;%# ... 
idx.bayesian = 0;%# % ...

%% parameters ... 
ps.audio_names = 'Kell_final_15_2.xlsx';
ps.conds = {'IT','TN','TT'};
ps.conds_new = {'S3','S2','S1'};
ps.downsample = 200; % 200; % original fs ... 
ps.frex4plot = 100;% maximum >>> ps.downsample/2;
step = 0.5; win = 2;
bins_start = 0:step:ps.frex4plot-win;
bins_end = bins_start+win;
ps.normalized_bins = [bins_start;bins_end];
ps.sig = 0.05;
ps.n_perm = 1000;
ps.n_points = 100; % how many data points to use in ttest ... 
ps.n_values = 30; % how many values were used to calculate each point ... 


function making_folder(folder_name)
if ~exist(folder_name,'dir')
    mkdir(folder_name);
end

