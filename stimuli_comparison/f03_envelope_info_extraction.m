function f03_envelope_info_extraction(paths,ps)
input_root_path = paths.stimuli;
output_root_path = paths.envelopes;

audio_files = dir(fullfile(input_root_path,'*.wav'));

%% extracting envelopes ...
n_str = fprintf('extracting envelopes for %5s (%3d of %3d) ...',0,0,0);
for audio_i = 1:length(audio_files)
    %% reading files ...
    tmp_audio_id = audio_files(audio_i).name;
    fprintf([repmat('\b',1,n_str),'extracting envelopes for %5s (%3d of %3d) ...'],tmp_audio_id(1:5),audio_i,length(audio_files));
    
    [tmp_audio_file,fs] = audioread(fullfile(input_root_path,tmp_audio_id));
    
    %% envelope extraction ...
    tmp_envelope = abs(hilbert(tmp_audio_file)); % sum(abs(y_phrase));
    tmp_envelope_ds = resample(tmp_envelope,ps.downsample,fs); % downsample ...
    
    %% plot signal and envelope ... 
%     h = figure;
%     subplot(311);
%     plot(tmp_audio_file);
%     xlim([0, length(tmp_audio_file)]);
%     xlabel('Time');
%     ylabel('Amplitude (a.u.)')
%     title('Speech signal');
%     
%     subplot(312);
%     plot(tmp_envelope_ds./max(tmp_envelope_ds),'r:', 'linew',1.5);
%     xlim([0, length(tmp_envelope_ds)]);
%     xlabel('Time');
%     ylabel('Amplitude (a.u.)')
%     title('Hilbert envelope (d.s. 200 Hz)');
%     
%     subplot(313);
%     addpath('E:\mpi\D\allToolBoxes\gammatone-master\auditory_toolkit');
%     gts = gammatonegram(tmp_audio_file,fs,0.025,0.010,512,1,fs/2,0);
%     plot(mean(gts,1)./max(mean(gts,1)),'k:','linew',1.5);
%     xlim([0, size(gts,2)]);
%     xlabel('Time');
%     ylabel('Amplitude (a.u.)')
%     title('Gammatone (512 bands)');
%     
%     close(h);
    %% fft envelope  ...
    signal_length = length(tmp_envelope_ds);
    f = ps.downsample.*(0:(signal_length/2))/signal_length;
    f=f(2:end);
    
    tmp_envelope_fft = 2*(abs(fft(tmp_envelope_ds)/length(tmp_envelope_ds)));
    tmp_envelope_fft = tmp_envelope_fft(2:length(f)+1);
    
    %% saving envelopes'frequency responses ...
    env_info(audio_i).name = tmp_audio_id(1:5);
    env_info(audio_i).audio_type = tmp_audio_id(1:2);
    env_info(audio_i).frequency = f';
    env_info(audio_i).envelope = tmp_envelope_ds;
    env_info(audio_i).envelope_fft = tmp_envelope_fft;
    
end
fprintf('\n');
%% reading audio name ...
all_audio_names = [];
for i=1:length(ps.conds)
    [~,~,tmp_cond_names]=xlsread(fullfile(input_root_path,ps.audio_names),ps.conds{i});
    tmp_cond_names(1,:)=[];
    all_audio_names = [all_audio_names;tmp_cond_names];
end
all_audio_names(:,1)=[];
[env_info.audio_name] = deal(all_audio_names{:});
save(fullfile(output_root_path,'envelope_info.mat'),'env_info');
