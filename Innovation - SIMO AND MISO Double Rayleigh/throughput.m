close all
clear all

N=10^4;
snr=-5:5:40;
EsN0=10.^(snr./10);
f=sqrt(0.5);
index=1;
sigma = 1; % the parameter
dSR = 1;
alpha = -4;
Ps = 1;

for k=-5:5:40
x=10^(k./10);
p=sqrt(1/x);
x1=randi([0,1]);
x2=2*x1-1;


h1=f*(randn(1,N) + j*randn(1,N));
h2=f*(randn(1,N) + j*randn(1,N));
g1= abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
g2= abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
h = 4.*((g1.*h1).*(h2.*g2));
n1=f*(randn(1,N) + j*randn(1,N));
n2=f*(randn(1,N) + j*randn(1,N));
n=n1.*n2;
y = h.*x2 + p.*n;
y1=h1.*x2 + p.*n1;
y2=h2.*x2 + p.*n2;


for kk=1:N
d(kk)=conj(h(kk)) * y(kk);
if(real(d(kk))>=0)
data_detect(kk)=1;
else
data_detect(kk)=0;

end


%For Variable Gain%
absolute_value = abs(h);
snr_th = 1.58489;
if(absolute_value(absolute_value > snr_th))
    absolute_value = absolute_value(absolute_value > snr_th);
else
end


%For Fixed Gain%
c1= 0.5*log2(1+EsN0);
if(c1(c1 < 2))
    out = c1(c1 < 2);
else
end

end


error = xor(x1,data_detect);
ber(index)=sum(error)/N;
snr(index)=k;
[snr(index) ber(index)];
N=N+1000;
index=index+1;
end

yf2 = [0.3679 0.6860 1.1087 1.4981 1.7063 1.8453 1.8810 1.9335 1.9826 2.0000];
yv2 = [0.2484 0.6951 1.4286 1.7536 1.9813 1.9845 1.9881 1.9952 1.9998 2.0000];
yf4 = [0.1671 0.4237 1.0121 1.8913 2.6854 3.2493 3.5441 3.8207 3.9607 3.9798];
yv4 = [0.0032 0.2183 0.9974 2.2370 3.1692 3.6267 3.8341 3.9450 3.9784 4.0000];
op1_sim2 = [0.3679 0.6860 1.1087 1.4981 1.7063 1.8453 1.8810 1.9335 1.9826 2.0000];
op2_sim2 = [0.2484 0.6951 1.4286 1.7536 1.9813 1.9845 1.9881 1.9952 1.9998 2.0000];
op1_sim4 = [0.1671 0.4237 1.0121 1.8913 2.6854 3.2493 3.5441 3.8207 3.9607 3.9798];
op2_sim4 = [0.0032 0.2183 0.9974 2.2370 3.1692 3.6267 3.8341 3.9450 3.9784 4.0000];

snr=-5:5:40;
plot(snr, yf2, '-b','Linewidth',2);  
hold on
plot(snr, yv2, '-kv','Linewidth',1);    
hold on
plot(snr, yf4, '-o','Linewidth',2);  
hold on
plot(snr, yv4, '-g','Linewidth',1);    
hold on
plot(snr, op1_sim2, 'mh','Linewidth',2);         
hold on
plot(snr, op2_sim2, '-r','Linewidth',2);
hold on
plot(snr, op1_sim4, '-y','Linewidth',2);         
hold on
plot(snr, op2_sim4, '-c','Linewidth',2);  
hold off

legend('Location','northwest')
LEG = legend('For R=2 Fixed Gain(Theory)', 'For R=2 Variable Gain(Theory)', 'For R=4 Fixed Gain(Theory)', 'For R=4 Variable Gain(Theory)', 'For R=2 Fixed Gain(Simulation)', 'For R=2 Variable Gain(Simulation)', 'For R=4 Fixed Gain(Simulation)','For R=4 Variable Gain(Simulation)');
LEG.FontSize = 8;
xlabel('Average SNR(dB)');
ylabel('Throughput');  

grid on
hold off
%quadgk
ylim([0.001 4]);
xlim([-5 40])
