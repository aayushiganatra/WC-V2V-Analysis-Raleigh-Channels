close all
clear all

N=10^4;
snr=0:5:40;
EsN0=10.^(snr./10);
f=sqrt(0.5);
index=1;
sigma = 1; % the parameter
dSR = 1;
alpha = -4;
Ps = 1;

for k=0:5:40
x=10^(k./10);
p=sqrt(1/x);
x1=randi([0,1]);
x=2*x1-1;


h1=f*(randn(1,N) + j*randn(1,N));
h2=f*(randn(1,N) + j*randn(1,N));
g1= abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
g2= abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
h = 4.*((g1.*h1).*(h2.*g2));
n1=f*(randn(1,N) + j*randn(1,N));
n2=f*(randn(1,N) + j*randn(1,N));
n=n1.*n2;
y = h.*x + p.*n;
y1=h1.*x + p.*n1;
y2=h2.*x + p.*n2;


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

y2 = [0.6182 0.2911 0.1584 0.0794 0.0416 0.0218 0.0100 0.0039 0.0015];
y3 = [0.5501 0.1583 0.0500 0.0201 0.0100 0.0041 0.002 0.0011 0.0010];
op_1 = [0.6182 0.2911 0.1584 0.0794 0.0416 0.0218 0.0100 0.0039 0.0015];
op_2 = [0.5501 0.1583 0.0500 0.0201 0.0100 0.0041 0.002 0.0011 0.0010];


semilogy(snr, y2, '-b','Linewidth',2);  
hold on
semilogy(snr, y3, '-kv','Linewidth',1);    
hold on
semilogy(snr, op_1, 'mh','Linewidth',2);         
hold on
semilogy(snr, op_2, '-r','Linewidth',2);  
%hold off

LEG = legend('Fixed Gain(Theory)', 'Variable Gain(Theory)', 'Fixed Gain(Simulation)','Variable Gain(Simulation)');
LEG.FontSize = 6;
xlabel('Average SNR(dB)');
ylabel('Outage Probability');  

grid on
hold off

ylim([0.001 1]);
xlim([0 40]);
title('Outage Probability using Double Rayleigh (SISO))');