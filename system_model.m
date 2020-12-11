close all
clear all

N=10^6; %number of bits or symbols
snr=0:5:40; %multiple Eb/No(SNR) value in db
f=sqrt(0.5);
EsN0=10.^(snr./10); % snr value(db) to linear scale
theory=0.5.*(1 - sqrt(EsN0./(EsN0+1)));
u = rand(N, 1); % generating uniform variates
sigma = 1; % the parameter
dSR = 1;
alpha = -4;
Ps = 1;
i=1;
c1= 0.5*log2(1+EsN0);
ok1 = dSR.^alpha;
ok = sqrt(ok1 * Ps);
out = 0;

for k=1:5:40
    out=0;
x=10^(k./10);
p=sqrt(1/x);
%mu=sqrt(x/(x+1));
%bera(k)=0.25.*(2-3.*mu + mu.^3); %Analytical formula for SIMO of 2 channels
x1=randi([0,1]); %Random generation of numbers
x2 = 2*x1-1;

%g1 = sigma * sqrt(-2 * log(u));
%g2 = sigma * sqrt(-2 * log(u));
h1=f*(randn(1,N) + j*randn(1,N));
h2=f*(randn(1,N) + j*randn(1,N));
%h = 4*(h1*g1 + h2*g2);
g1= abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
g2= abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
h = 4.*((g1.*h1).*(h2.*g2));
n1=f*(randn(1,N) + j*randn(1,N));
n2=f*(randn(1,N) + j*randn(1,N));
n=n1.*n2;
y1=h.*x2 + p.*n;


for kk=1:N
b(kk)=conj(h(kk)) * y1(kk); %Calculation for real and imaginary parts of signals
if(real(b(kk))>=0) %inphase demodulation
data_detect_b(kk)=1; %Detection of real part for BER
else
data_detect_b(kk)=0;
end
absolute_value = abs(h);
snr_th = 1.58489;
if(c1 > 1)
  out = out + 1;
   else
end
end

sout(i)=out./N;
i=i+1;

%error_b = xor(x1,data_detect_b); %Bit error rate for 1*1 SISO wireless system
%bers_b(k)=sum(error_b)/N; %Sum of errors by the total transmission bits

    %for (snr=0:5:40);
    %op=1-(exp(-snr./bers_b)); 
     %end   
end


y2 = [0.55 0.48 0.34 0.1 0.078 0.037 0.01 0.006 0.002];
y3 = [0.48 0.33 0.089 0.045 0.01 0.004 0.002 0.0012 0.000];
op_1 = [0.55 0.48 0.34 0.1 0.078 0.037 0.01 0.006 0.002];
op_2 = [0.48 0.33 0.089 0.045 0.01 0.004 0.002 0.0012 0.000];

semilogy(snr,sout(1:26),'g*')
%semilogy(snr, y2, '-b','Linewidth',2);  
%hold on
%semilogy(snr, y3, '-kv','Linewidth',1);    
%hold on
%semilogy(snr, op_1, 'mh','Linewidth',2);         
%hold on
%semilogy(snr, op_2, '-r','Linewidth',2);  
%hold off

LEG = legend('Fixed Gain(Theory)', 'Variable Gain(Theory)', 'Fixed Gain(Simulation)','Variable Gain(Simulation)');
LEG.FontSize = 6;
xlabel('Average SNR(dB)');
ylabel('Outage Probability');  

grid on
hold off

ylim([0.001 1]);
xlim([0 40])