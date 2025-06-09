function f02_stimuli_spectrogram_plot(paths,ps)
%%
input_root_path = paths.stimuli;
output_root_path = paths.audio_spectrogram;

% [~,~,raw] = xlsread(fullfile(input_root_path,ps.audio_names),'IT');
audiofiles = dir(fullfile(input_root_path,'*.wav'));

for cond_i = 1:length(ps.conds)
    tmp_cond = ps.conds{cond_i};
    tmp_cond_idx = cellfun(@(x) ~isempty(x),regexpi({audiofiles.name},[tmp_cond,'\w*']),'UniformOutput',1);
    tmp_conds_audiofiles = audiofiles(tmp_cond_idx);
    feature_plot(tmp_cond,tmp_conds_audiofiles,input_root_path,output_root_path,ps);
end


%% actual plot ...
function feature_plot(cond,audiofiles_info,input_path,output_path,ps)
%% reading audio names ...
[~,~,audio_names] = xlsread(fullfile(input_path,ps.audio_names),cond);
audio_names(1,:)=[];
for audio_i = 1:length(audiofiles_info)
    %% reading files ...
    tmp_audio_name = [strrep(audiofiles_info(audio_i).name(1:5),'_',''),': ',audio_names{audio_i,2}];
    [tmp_audio,fs] = audioread(fullfile(input_path,audiofiles_info(audio_i).name));
    
    %%
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    %     h = figure('Position',[100 100 800 800]);
    
    
    %% plot spectrogram ...
    subplot(122)
    spectrogram(tmp_audio,256,250,[],fs,'yaxis');
    set(gca,'XTick',[0:0.5:2]);
    set(gca,'XTickLabel',strsplit(num2str([0:500:2000])));
    %     colorbar('off');
    colormap('jet');
    axis square; set(gca,'clim',[-130 -30]);
    title(tmp_audio_name); xlabel('Time (ms)')
    
    
    %% plot envelope ...
    subplot(221)
    tmp_envelope = abs(hilbert(tmp_audio)); % sum(abs(y_phrase));
    tmp_envelope = resample(tmp_envelope,ps.downsample,fs); % downsample ...
    
    p1 = plot([0:length(tmp_envelope)-1].*(1/ps.downsample),tmp_envelope,'linew',1.2);
    xlim([0,length(tmp_envelope)-1].*(1/ps.downsample)); ylim([0 0.5])
    set(gca,'XTick',[0:0.5:2]);
    set(gca,'XTickLabel',strsplit(num2str([0:500:2000])));
    xlabel('Time (ms)');
    ylabel('Amplitude (a.u.)');
    legend(p1,['Envelope: ',strrep(audiofiles_info(audio_i).name(1:5),'_','')]);
    title(['temporal envelopes ', '(D.S. ',sprintf('%.3d',ps.downsample),' Hz)']);
    axis square;
    
    %% plot fft envelope  ...
    signal_length = length(tmp_envelope);
    f = ps.downsample.*(0:(signal_length/2))/signal_length;
    f=f(2:end);
    
    subplot(223)
    tmp_audio_amp = 2*(abs(fft(tmp_envelope)/length(tmp_envelope)));
    tmp_audio_amp = tmp_audio_amp(2:length(f)+1);
    
    %     tmp_baseline = max(tmp_audio_amp);
    %     tmp_audio_amp = 20*log10(tmp_audio_amp./tmp_baseline);
    plot_idx = find(f<=ps.frex4plot);
    plot(f(plot_idx),tmp_audio_amp(plot_idx),'linew',1.0); ylim([0 0.1]);
    legend(strrep(audiofiles_info(audio_i).name(1:5),'_',''));
    axis square; box on;
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (a.u.)');
    title('response spectrum of envelope');
    
    %% saving pics ...
    tmp_pic_name = strrep([strrep(audiofiles_info(audio_i).name(1:5),'_',''),' ',audio_names{audio_i,2}],' ','_');
    tmp_pic_name = tmp_pic_name(1:find(ismember(tmp_pic_name,'.')==1)-1);
    pic_saving_path = fullfile(output_path,cond);
    if ~exist(pic_saving_path,'dir')
        mkdir(pic_saving_path);
    end
    print(gcf,fullfile(pic_saving_path,tmp_pic_name),'-dtiff');
    close(h);
end
