close all;

SNR = -5:5:35;
R2 = 2;
R4 = 4;

yf = [0.95 0.93 0.78 0.13 0.078 0.037 0.012 0.006 0.002];
yv = [0.98 0.88 0.289 0.045 0.01 0.004 0.002 0.0012 0.000];

Tf2 = R2*(1 - yf)
Tv2 = R2*(1 - yv)

Tf4 = R4*(1 - yf)
Tv4 = R4*(1 - yv)

plot(SNR,Tf2,'-b','Linewidth', 2)
hold on
plot(SNR,Tv2,'-r','Linewidth', 2)
hold on 
plot(SNR,Tf4,'-o','Linewidth', 2)
hold on
plot(SNR,Tv4,'-g','Linewidth', 2)
legend( 'Fixed Gain(Theory) for R=2', 'Variable Gain(Theory) for R=2', 'Fixed Gain(Theory) for R=4', 'Variable Gain(Theory) for R=4');
xlabel('Average SNR(dB)');
ylabel('Throughput');

grid on
hold off
ylim([0 5]);
xlim([-5 40]);