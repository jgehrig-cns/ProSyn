function threshold = f05_temporal_similarity_on_envelopes(paths,ps)
rng('shuffle')
input_root_path = paths.envelopes;
output_root_path = input_root_path;
load(fullfile(input_root_path,'envelope_info.mat'));

%% cross-correlation ...
max_cross_corr = zeros(length(env_info),length(env_info));
all_envelopes = {env_info.envelope}';
for i=1:length(all_envelopes)
    tmp_env1 = all_envelopes{i};
    for j=1:length(all_envelopes)
        tmp_env2 = all_envelopes{j};
        max_cross_corr(i,j) = max(zscore(xcorr(tmp_env1,tmp_env2)));
    end
end

%% permutation test ...
for i=1:ps.n_perm
    tmp_env = all_envelopes{randperm(length(all_envelopes),1)};
    tmp_rand_env = Shuffle(all_envelopes{randperm(length(all_envelopes),1)});
    tmp_rand_env = rand(length(tmp_rand_env),1).*tmp_rand_env; % this operation is applausible if you did stimuli's amplitude normalization ...
    null_corr(i) = max(zscore(xcorr(tmp_env,tmp_rand_env)));
end
null_corr = sort(null_corr)';
threshold = null_corr(ps.n_perm*(1-ps.sig));

%% plot xcorr results ...
max_cross_corr_padded = max_cross_corr;
max_cross_corr_padded(length(max_cross_corr_padded)+2,length(max_cross_corr_padded)+2)=0;

% rearrange columns ... 
new_col_order = [1:40,121,41:80,122,81:120];
max_cross_corr_padded = max_cross_corr_padded(:,new_col_order);

% rearrange rows ... 
new_row_order = new_col_order;
max_cross_corr_padded = max_cross_corr_padded(new_row_order,:);
save(fullfile(output_root_path,'similarity_data.mat'),'max_cross_corr_padded','threshold'); % save it for statistic analysis ... 

max_cross_corr_padded(max_cross_corr_padded<=threshold)=0;

% apply transparency ... 
alpha_degree = 0.1;
trans_mat = ones(length(max_cross_corr),length(max_cross_corr));
trans_mat(121:122,:)=alpha_degree;
trans_mat(:,121:122) = alpha_degree;
trans_mat = trans_mat(:,new_col_order);
trans_mat = trans_mat(new_row_order,:);

%% plot ... 
h = figure;
clim_range = [0 3];
imagesc(max_cross_corr_padded,'AlphaData',trans_mat);
set(gca,'ydir','normal',...
    'clim',clim_range,...
    'XTick',[20,61,102],...
    'XTickLabel',ps.conds_new,...
    'YTick',[20,61,102],...
    'YTickLabel',ps.conds_new,...
    'FontWeight','bold');
colormap('jet');axis square;
title('RSM','FontWeight','bold');
% cc =colorbar('Ticks',[clim_range(1), threshold, clim_range(2)],...
%     'TickLabels',{num2str(clim_range(1)),['\theta: ',sprintf('%.2f',threshold)],num2str(clim_range(2))});
cc= colorbar;
cc.Label.String = 'Similarity (z-score)';

%% saving the picture ... 
print(gcf,fullfile(output_root_path,'similarity_comparison_on_envelopes'),'-dtiff');
close(h);
