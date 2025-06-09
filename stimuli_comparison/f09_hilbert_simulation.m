close all;
clear;clc;

%%  
t = (1:131072)*(1e-4);
gt = sin(t).*sin(18*t);
envelope = abs(hilbert(gt));
figure; 
plot(t,gt,'k-','linew',1.5);
hold on;
plot(t,envelope,'r:','linew',1.5);

plot(t,imag(hilbert(gt)),'b--','linew',1.5);
xlim([0.25,12.5]);
xlabel('Time');
ylabel('Amplitude (a.u.)');
title('Hilbert Tranform');
legend('Signal','Envelope','Hilbert');
