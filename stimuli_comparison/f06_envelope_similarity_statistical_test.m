function f06_envelope_similarity_statistical_test(paths,ps,threshold)
rng('shuffle');


input_root = paths.envelopes;
output_root = input_root;

load(fullfile(input_root,'similarity_data'));

it_vs_tt = max_cross_corr_padded(1:40,83:122); % ## c1 ...
% it_vs_tt(it_vs_tt<threshold)=0;
% imagesc(tn_vs_tt),set(gca,'YDir','normal'),colormap('jet'),axis square; % checking correctness ...

it9tn = max_cross_corr_padded(1:40,42:81); % ## c2 ...
% it_vs_tn(it_vs_tn<threshold)=0;

tn_vs_tt = max_cross_corr_padded(41:80,83:122); % ## c3 ...
% tn_vs_tt(tn_vs_tt<threshold)=0;

[a,b,c,d] = ttest(it_vs_tt(:), tn_vs_tt(:));

%% statistics generation ...
data4stats = zeros(ps.n_points,3);
for i=1:ps.n_points
    tmp_it_vs_tt = mean(it_vs_tt(randperm(length(it_vs_tt(:)),ps.n_values)));
    tmp_it_vs_tn = mean(it9tn(randperm(length(it9tn(:)),ps.n_values)));
    tmp_tn_vs_tt = mean(tn_vs_tt(randperm(length(tn_vs_tt(:)),ps.n_values)));
    data4stats(i,1:3) = [tmp_it_vs_tt,tmp_it_vs_tn,tmp_tn_vs_tt];
end

%% 1-way anova ...
[p_anova,~,~] = anova1(data4stats,[],'off');
% [p_anova,~,~] = anova1(data4stats); % for writing ... 


%% independent comparisons ...
[h_c1_vs_c2,p_c1_vs_c2,~,stats_c1_vs_c2] = ttest(data4stats(:,1),data4stats(:,2));
[h_c2_vs_c3,p_c2_vs_c3,~,stats_c2_vs_c3] = ttest(data4stats(:,2),data4stats(:,3));
[h_c1_vs_c3,p_c1_vs_c3,~,stats_c1_vs_c3] = ttest(data4stats(:,1),data4stats(:,3));

%% print results ...
if p_anova<0.05
    fprintf('c1 (it&tt) V.S. c2 (it&tn), ( h = %d, p = %s, df = %d , stats = %.4f)... \n\n',h_c1_vs_c2 ,sprintf('%0.3e',p_c1_vs_c2/3),stats_c1_vs_c2.df,stats_c1_vs_c2.tstat); % /3 >> Bonf correction ...
    fprintf('c2 (it&tn) V.S. c3 (tn&tt), ( h = %d, p = %s, df = %d , stats = %.4f)... \n\n',h_c2_vs_c3 ,sprintf('%0.3e',p_c2_vs_c3/3),stats_c2_vs_c3.df,stats_c2_vs_c3.tstat);
    fprintf('c1 (it&tt) V.S. c3 (tn&tt), ( h = %d, p = %s, df = %d , stats = %.4f)... \n\n',h_c1_vs_c3 ,sprintf('%0.3e',p_c1_vs_c3/3),stats_c1_vs_c3.df,stats_c1_vs_c3.tstat);
else
    fprintf('failed to find any significant difference ... \n\n');
end

%% plotting ...
h1 = figure;
colors = [0.35 0.35 0.35;
    0.6 0.6 0.6;
    0.9 0.9 0.9];
for i=1:3
    p(i) = bar(i,mean(data4stats(:,i)),'FaceColor',colors(i,:),'EdgeColor',colors(i,:),'LineWidth',1.5,'BarWidth',0.6);
    se(i) = std(data4stats(:,i))/sqrt(length(data4stats)-1);
    hold on;
    plot([i i],[mean(data4stats(:,i))-se(i), mean(data4stats(:,i))+se(i)],'color',[0 0 0],'linew',2);
    sig_pos(i) = mean(data4stats(:,i))+se(i);
end
text_x = {'IT&TT','IT&TN','TN&TT'};
% legend('boxoff');
ylabel('Similarity (a.u.)','fontsize',10,'fontweight','bold');
set(gca,'XTick',1:3,'XTickLabel',text_x,'YTick',2:0.04:2.2);
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
title('Similarity comparison across conditions');
axis square; box on;

%% adding significance ...
line1_y = sig_pos(1)+0.01;
plot(1:2,repmat(line1_y,1,2),'color','k','linestyle',':','linew',1.5);
text(1.5,line1_y+0.0025,'***','HorizontalAlignmen','center','fontsize',10,'fontweight','bold');

line2_y = sig_pos(2)+0.01;
plot(2:3,repmat(line2_y,1,2),'color','k','linestyle',':','linew',1.5);
text(2.5,line2_y+0.0025,'***','HorizontalAlignmen','center','fontsize',10,'fontweight','bold');

line3_y = line1_y  + (line1_y-line2_y);
plot(1:3,repmat(line3_y,1,3),'color','k','linestyle',':','linew',1.5);
text(2,line3_y+0.0025,'***','HorizontalAlignmen','center','fontsize',10,'fontweight','bold');

ylim([2 line3_y+0.05]);
%% saving the figure ...
legend(p,text_x);

print(fullfile(output_root,'envelope_temporal_comparison'),'-dtiff','-r600');
close(h1);

%% additional comparison 2023/08/11 ...
tt9it = data4stats(:,1);
[h,p,~,stats] = ttest(tt9it,threshold);
if p/3<0.05
    h=1;
end
fprintf('tt&it against theta, ( h = %d, p = %s, df = %d , stats = %.4f)... \n\n',...
    h ,sprintf('%0.3e',p/3),stats.df,stats.tstat); % /3 >> Bonf correction ...


% 1000 times ttest ... 
for i=1:10000
    tmp_tt9tn = [];
    for j=1:ps.n_points
        k = mean(tn_vs_tt(randperm(length(tn_vs_tt(:)),ps.n_values)));
        tmp_tt9tn(j) = k;
    end
%     [h,p,~,stats] = ttest(tmp_tt9tn,threshold*ones(size(tmp_tt9tn)),'tail','right');
    [h,p,~,stats] = ttest(tmp_tt9tn,threshold,'tail','right');
    a(i,1:2) = [p/3,stats.tstat];
end
p = mean(a(:,1));
if p<0.05
    h=1;
else
    h=0;
end
fprintf('tt&tn against theta, ( h = %d, p = %s, df = %d , stats = %.4f)... \n\n',...
    h ,sprintf('%0.3e',p),stats.df,mean(a(:,2))); % /3 >> Bonf correction ...
a=1;



