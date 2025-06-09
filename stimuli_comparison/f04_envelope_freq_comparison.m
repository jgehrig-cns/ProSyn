function f04_envelope_freq_comparison(paths,ps)
input_root_path = paths.envelopes;
load(fullfile(input_root_path,'envelope_info.mat'));

%% extracting normalized frequency response ...
norm_frex = mean(ps.normalized_bins,1);
env_compare_matrix = [];
n_str = fprintf('normalize envelope frequency bin %3d of %3d ...',0,0);
for env_i=1:length(env_info)
    env_info(env_i).norm_freq = norm_frex;
    fprintf([repmat('\b',1,n_str),'normalize envelope frequency bin %3d of %3d ...'],env_i,length(env_info));
    for bin_i=1:length(norm_frex)
        tmp_freq_range = ps.normalized_bins(:,bin_i);
        tmp_f_idx = [env_info(env_i).frequency]>=tmp_freq_range(1)&[env_info(env_i).frequency]<=tmp_freq_range(2);
        tmp_amp(bin_i) = mean( env_info(env_i).envelope_fft(tmp_f_idx));
    end
    env_info(env_i).norm_amp = tmp_amp';
end
fprintf('\n');

%% plot normalized freqquency response ...
color4env = {'r','g','b'};
% color4env = {[0 0 0],[221 0 0]./255,[255 206 0]./255};

h=figure;
for cond_i = 1:length(ps.conds)
    tmp_cond_idx = strcmpi({env_info.audio_type},ps.conds{cond_i});
    tmp_cond_env_info = env_info(tmp_cond_idx);
    data4statistics(:,:,cond_i) = cat(2,tmp_cond_env_info.norm_amp);
    [frex4plot,tmp_avg_env,tmp_shading_x,tmp_shading_env] = plot_info(tmp_cond_env_info);
    p(cond_i) = plot(frex4plot,tmp_avg_env,'color',color4env{cond_i},'linew',1.2);
    hold on;
    fill(tmp_shading_x,tmp_shading_env,color4env{cond_i},'EdgeColor','none','FaceAlpha',0.2);
end
legend(p,ps.conds);
title('Response spectrum of envelopes');
xlabel('Frequency (Hz)','fontsize',10,'fontweight','bold');
ylabel('Amplitude (a.u.)','fontsize',10,'fontweight','bold');
axis square;
ax= gca; 
ax.YAxis.Exponent = -2;

%% adding significance ...
ylim_range = [-inf 0.05];
% sig_frex = [2,24.5,25,91.5,92,92.5,93,98.5];
sig_frex = [2];

for i=1:length(sig_frex)
    s = plot([sig_frex(i) sig_frex(i)], [0.02 0.045],'linestyle',':','linew',1.2,'color',[0 0 0 0.5]);
    hold on;
end
xlim([1 inf]);
ylim(ylim_range);
text(2.2,0.043,'f = 2 Hz','color','k',...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle',...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontAngle','italic');

set(gca,'XScale','linear');

%% saving the picture ... 
print(gcf,fullfile(input_root_path,'envelope_spectrum'),'-dtiff','-r600');
close(h);

%% statistics ...
for i=1:size(data4statistics,1)
    tmp_data = squeeze(data4statistics(i,:,:));
    conds = repmat(ps.conds,size(tmp_data,1),1);
    p_info(i) = anova1(tmp_data(:),conds(:),'off');
end
[ind, thres] = FDR(p_info, ps.sig); % FER correction ...

if isempty(ind)
    fprintf('one-way ANOVA failed to find any difference across conditions ... \n');
end

%% saving data for Bayesian statistics ... 
save(fullfile(input_root_path,'info4bayesian.mat'),'norm_frex','p_info','data4statistics');

%% plot statistics ... 
h1 = figure;
scatter(norm_frex,p_info,'filled');
hold on;
plot([0 100],0.05*ones(size([0 100])),'r:','linew',2);
text(50,0.1,'p = 0.05','color','k');
hold off; box on; axis square; 
xlabel('Frequency (Hz)');
ylabel('P-value');
title('one-way ANOVA (un-corrected)');
print(gcf,fullfile(input_root_path,'env_freq_statistics'),'-dtiff','-r600');
close(h1);



%%
function [frex4plot,avg_env,shading_x,shading_env] = plot_info (env_data)
frex = env_data(1).norm_freq;
envelopes = cat(2,env_data.norm_amp); 
freq_cutoff = 10; % in Hz ... 
freq_cutoff_idx = dsearchn(frex',freq_cutoff);

frex4plot = frex(1:freq_cutoff_idx);
env4plot = envelopes(1:freq_cutoff_idx,:);

avg_env = mean(env4plot,2)';
se_env = 2*std(env4plot,0,2)'./(sqrt(length(env_data)-1));


shading_x = [frex4plot,fliplr(frex4plot)];
shading_env = [avg_env+se_env,fliplr(avg_env-se_env)];

%% FDR correction ...
function [ ind, thres ] = FDR( p_list, alpha, corrected )
% Computes the False Discovery Rate according to Benjamini and Hochberg (1995).
%
% Inputs:
% p_list - list of p values
% alpha - the desired alpha threshold. Default: 0.05
% corrected - set to true if correction for dependencies is to be applied, according to Benjamini
% and Yekutieli (2001) (this is probably not the common case).
%
% outputs:
% ind - the indexes of significant p-values within p_list
% thres - the p-value which served as the actual threshold in this test.
%
% Written by Edden Gerber, lab of Leon Y. Deouell, 2012
% Please send bug reports and requsts to edden.gerber@gmail.com
%
n_vals = length(p_list);
num_tests = n_vals; % there was some reason that in some cases you may want to set this to
% a lower value, but I don't remember what it was.
if nargin < 2
    alpha = 0.05;
end
if nargin < 3
    corrected = false;
end
p_sorted = sort(p_list,'descend');
if corrected
    comp = (num_tests:-1:1)/num_tests * alpha / sum((1:num_tests)/num_tests);
else
    comp = (num_tests:-1:1)/num_tests * alpha;
end
comp = comp((end-n_vals+1):end);
i = find(p_sorted <= comp,1,'first');
if isempty(i)
    thres = 0;
else
    thres = p_sorted(i);
end
ind = find(p_list<=thres);

