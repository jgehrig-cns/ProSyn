%% 
close all;
clear;clc;

%% root_path ...
root_path = 'E:\mpi\K\projects_mpi\p03_frankfurt_EcoG_Done';

%% paths and idx setting ...
[paths,idx, ps] = f01_path_and_idx_setting(root_path);

%% plot spectrogram ...
if idx.plot_spectrogram
    f02_stimuli_spectrogram_plot(paths,ps);
end

%% envelope extraction ...
if idx.envelope_info_extraction
    f03_envelope_info_extraction(paths,ps);
end

%% frequency domain comparison ...
if idx.envelope_freq_statistics
    f04_envelope_freq_comparison(paths,ps);
end

%% time domain comparison ...
if idx.RSM
    threshold = f05_temporal_similarity_on_envelopes(paths,ps);
end

%% time domain significance test ...
if idx.time_domain_statistics
    f06_envelope_similarity_statistical_test(paths,ps,threshold);
end

%% bayesian inference on spectrum of envelope ...
if idx.bayesian
    f07_env_fft_bayesian(paths,ps);
end