close all
clear all

N=100000;
snr=0:5:40; %multiple Eb/No(SNR) value in db
f=sqrt(0.5);
index=1;
M = 2;
g = sin(pi/M)^2;
sigm = 1; % channel variance
sigma = 1;

Mr=1;
Mt=100000;


for k=1:5:40
x=10^(k./10);
p=sqrt(1/x);
x1=randi([0,1]);
x=2*x1-1;


h1=f*(randn(1,N) + j*randn(1,N));
h2=f*(randn(1,N) + j*randn(1,N));
%h = 4*(h1*g1 + h2*g2);
g1= abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
g2= abs(sigma*randn(1,N)+1i*sigma*randn(1,N));
h = 4.*((g1.*h1).*(h2.*g2));

n1=f*(randn(1,N) + j*randn(1,N));
n2=f*(randn(1,N) + j*randn(1,N));
n=n1.*n2;
y1=h.*x + p.*n;

for kk=1:N
   
     b(kk)=conj(h(kk)) * y1(kk);           %Calculation for real and imaginary parts of signals
  if(real(b(kk))>=0)                        %inphase demodulation
     data_detect_b(kk)=1;                   %Detection of real part for BER
 else
     data_detect_b(kk)=0;
  end
end
 %For Fixed Gain%
    c1= 0.5*log2(1+x);

   if (c1 < 3)
       out = c1;
   end    
error_b = xor(x1,data_detect_b);    %Bit error rate for 1*1 SISO wireless system
bers_b(k)=sum(error_b)/N;          %Sum of errors by the total transmission bits

end

y5=[0.00980582617584078 0.00369323661384938 0.001243758592258976 0.0004178433253485600 0.00013771151735 0.000085985548944 0.00003144982824673 0.000018402393774 0.00000513194026 ]
op_5 = [0.00980582617584078 0.00369323661384938 0.001243758592258976 0.0004178433253485600 0.00013771151735 0.000085985548944 0.00003144982824673 0.000018402393774 0.00000513194026 ]


semilogy(snr, y5, '-r','Linewidth',2);
hold on
semilogy(snr, op_5, 'mh','Linewidth',2);         
hold off

title('Outage Probability Using Double Rayleigh(SIMO)')
LEG = legend('SIMO(Theory)','SIMO(Simulation)');
LEG.FontSize = 5;
xlabel('Average SNR(dB)');
ylabel('Outage Probability');
grid on
hold off
ylim([0.0000001 1]);
xlim([0 40])
