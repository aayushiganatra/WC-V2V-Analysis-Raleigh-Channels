close all
clear all

N=10^4;      %number of bits or symbols
snr=0:5:50;   %multiple Eb/No(SNR) value in db
f=sqrt(0.5);  
EsN0=10.^(snr./10);    % snr value(db) to linear scale
theory=0.5.*(1 - sqrt(EsN0./(EsN0+1)));
%u = rand(N, 1); % generating uniform variates
sigma = 1; % the parameter
dSR = 0.5;
alpha = -4;
Ps = 1;


for k=1:1:50
x=10^(k./10);          
p=sqrt(1/x);            
%mu=sqrt(x/(x+1));
%bera(k)=0.25.*(2-3.*mu + mu.^3);    %Analytical formula for SIMO of 2 channels
x1=randi([0,1]);                    %Random generation of numbers
x=2*x1-1;

ok1 = dSR.^alpha;
ok = sqrt(ok1 * Ps);
g1 = abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
g2 = abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
%g1 = sigma * sqrt(-2 * log(u));
%g2 = sigma * sqrt(-2 * log(u));
%h1=f*(randn(1,N) + j*randn(1,N));
%h2=f*(randn(1,N) + j*randn(1,N));
h = ok.*(g1.*g2);
n1=f*(randn(1,N) + j*randn(1,N));  
n2=f*(randn(1,N) + j*randn(1,N));
n=n1.*n2
y1=h.*x + p.*n;  

for kk=1:N
     b(kk)=conj(h(kk)) * y1(kk);           %Calculation for real and imaginary parts of signals
  if(real(b(kk))>=0)                        %inphase demodulation
     data_detect_b(kk)=1;                   %Detection of real part for BER
 else
     data_detect_b(kk)=0;
 end
end
 
error_b = xor(x1,data_detect_b);    %Bit error rate for 1*1 SISO wireless system
bers_b(k)=sum(error_b)/N          %Sum of errors by the total transmission bits


end
y2=[0.316 0.2511 0.145 0.100 0.06 0.02511 0.0125 0.00594 0.00316 0.002 0.00158];
y3=[0.316 0.2000 0.080 0.050 0.02 0.00711 0.005  0.00300 0.0015 0.0010 0.0009];

snr = linspace(0, 50);
y2 = linspace(0.316, 0.00158);
y3 = linspace(0.316, 0.0009);
bers_b = linspace(1,N);

semilogy(snr, y2, '-o','Linewidth',2);  
hold on
semilogy(snr, y3, '-g','Linewidth',2);    
hold on
%semilogy(snr, theory, '-r','Linewidth',2);         %Plotting 1*1 Analytical
%hold on
%'-.*g'
semilogy(snr, bers_b, '-b', 'Linewidth',1);       %Plotting 1*1 Simulation
hold on
legend( '1*1 analytical', '1*1 simulation', '1*2 analytical', '1*2 simulation', '1*4 simulation');
title(' (1*2) vs (1*4) SIMO wireless system');
xlabel('SNR(db)');
ylabel('SER');  

grid on
%hold off
ylim([0.0001 1]);
xlim([0 50]);