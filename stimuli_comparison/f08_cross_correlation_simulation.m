close all;
clear; clc; 
%% generate sequences ... 
t = 0:0.05:1;
n = length(t);

y1 = [zeros(1,randperm(n,1)),sin(2*pi*t),zeros(1,randperm(n,1))];
y2 = [zeros(1,randperm(n,1)),sin(2*pi*t),zeros(1,randperm(n,1))];
y3 = [zeros(1,randperm(n,1)),sin(2*1.7*pi*t),zeros(1,randperm(n,1))];

y4 = rand(1,length(y1));
y5 = rand(1,length(y2));
sm = xcorr(y1,y2);

%% plot sine waves ... 
h=figure; 
subplot(221)
plot(1:length(y1),y1,'k','linew',1.5);
hold on;
plot(1:length(y2),y2, 'r:','linew',1.5);
plot(1:length(y3),y3, 'b--','linew',1.5)


legend('y1','y2','y3');
title('Sine waves');
xlabel('Time');
ylabel('Amplitue (a.u.)');
hold off;

subplot(223)
[a,~] = xcorr(y1,y2);
plot(zscore(a),'r:','linew',1.0)
hold on;
[b,~] = xcorr(y1,y3);
plot(zscore(b),'b--','linew',1.0)
legend('y1&y2','y1&y3','Location','best')

xlabel('Lag')
ylabel('Similarity (z-score)');
title('Cross correlation');
ylim([-2,4]);

%% plot random signals ... 
subplot(222)
plot(1:length(y4),y4,'k','linew',1.5);
hold on;
plot(1:length(y5),y5, 'r:','linew',1.5);
legend('y1-rand','y2-rand');
title('Random signals');
xlabel('Time');
ylabel('Amplitue (a.u.)');
hold off;

subplot(224)
[a,b] = xcorr(y4,y5);
plot(zscore(a),'r:','linew',0.8)
legend('y1&y2','Location','best')

xlabel('Lag')
ylabel('Similarity (z-score)');
title('Cross correlation');
ylim([-2,4]);
