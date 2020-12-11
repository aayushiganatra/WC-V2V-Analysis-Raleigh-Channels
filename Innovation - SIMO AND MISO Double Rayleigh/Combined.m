close all
clear all

N=100000;                                               %Total number of bits or symbols
snr = 0:5:40;                                          %Multiple Es/N0 values that are in dB
f=sqrt(0.5);
Mr=1;
Mt=100000; 
sigma = 1;

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


for ii = 1:length(snr)
    %From the Transmitter side
    x_eq = rand(1,N) > 0.5;                               %Generation of 0,1 values with equal probability
    x1 = 2*x_eq-1;                                         %BPSK modulation
   
    %Creating the Alamouti STBC
    input_arg = zeros(2,N);
    input_arg(:,1:2:end) = (1/sqrt(2))*reshape(x1,2,N/2);    
    input_arg(:,2:2:end) = (1/sqrt(2))*(kron(ones(1,N/2),[-1;1]).*flipud(reshape(conj(x1),2,N/2)));
    channelCoeff = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)];                         %Channel coefficient of Rayleigh channel
    channelCoefficientMod = kron(reshape(channelCoeff,2,N/2),ones(1,2));              
    noise = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)];                                      %White gaussian noise
   
    receiver = sum(channelCoefficientMod.*input_arg,1) + 10^(-snr(ii)/20)*noise;
    %From the Receiver side
    receiverMod = kron(reshape(receiver,2,N/2),ones(1,2));                              
    receiverMod(2,:) = conj(receiverMod(2,:));                                          
 
    %Forming the Equalization matrix
    equalizationMatrix_eq = zeros(2,N);
    equalizationMatrix_eq(:,[1:2:end]) = reshape(channelCoeff,2,N/2);                
    equalizationMatrix_eq(:,[2:2:end]) = kron(ones(1,N/2),[1;-1]).*flipud(reshape(channelCoeff,2,N/2));
    equalizationMatrix_eq(1,:) = conj(equalizationMatrix_eq(1,:));                            
    eqMatrixPower = sum(equalizationMatrix_eq.*conj(equalizationMatrix_eq),1);
    receiverHat_1 = sum(equalizationMatrix_eq.*receiverMod,1)./eqMatrixPower;                
    receiverHat_1(2:2:end) = conj(receiverHat_1(2:2:end));
   
    %Taking real part from the Receiver Side
    xHat = real(receiverHat_1)>0;
   
    %Total counting of the errors
    totalErrors(ii) = size(find([x_eq- xHat]),2);
  
    
    %For Fixed Gain%
    c1= 0.5*log2(1+x);

   if (c1 < 3)
       out = c1;
   end    
end
end

simulationAlamouti = totalErrors/N;                                                     %Simulation of BER
EbN0Linear = 10.^(snr/10);

%SISO Double Rayleigh
y1 = [0.6182 0.2911 0.1584 0.0794 0.0416 0.0218 0.0100 0.0039 0.0015];
op_1 = [0.6182 0.2911 0.1584 0.0794 0.0416 0.0218 0.0100 0.0039 0.0015];


%MISO Double Rayleigh
y2=4*[0.054921000000000  0.0239440000000000 0.011580990000000  0.004298000000000000  0.0020141000000000000 0.0010000000000000  0.00050000000000 0.0001811 0.00012511 ];
op_2 = 4*[0.054921000000000  0.0239440000000000 0.011580990000000  0.004298000000000000  0.0020141000000000000 0.0010000000000000  0.00050000000000 0.0001811 0.00012511 ];

%SIMO Double Rayleigh
y3=6*[0.00980582617584078 0.00369323661384938 0.001243758592258976 0.0004178433253485600 0.00013771151735 0.000085985548944 0.00003144982824673 0.000018402393774 0.00000513194026 ];
op_3 = 6*[0.00980582617584078 0.00369323661384938 0.001243758592258976 0.0004178433253485600 0.00013771151735 0.000085985548944 0.00003144982824673 0.000018402393774 0.00000513194026 ];

%SISO Simulation and Analytical
semilogy(snr, y1, '-b','Linewidth',2);  
hold on
semilogy(snr, op_1, 'mh','Linewidth',2);         
hold on
%MISO Simulation and Analytical
semilogy(snr, y2, '-b','Linewidth',2);  
hold on
semilogy(snr, op_2, '-kv','Linewidth',1);    
hold on
%SIMO Simulation and Analytical
semilogy(snr, y3, '-r','Linewidth',2);
hold on
semilogy(snr, op_3, 'mh','Linewidth',2);         
hold off

title('Outage Probability Using Double Rayleigh')
LEG = legend('Fixed Gain(Theory)', 'Fixed Gain(Simulation)', 'MISO(Fixed Gain)','MISO(Simulation)','SIMO(Fixed Gain)','SIMO(Simulation)');
LEG.FontSize = 5;
xlabel('Average SNR(dB)');
ylabel('Outage Probability');
grid on
hold off
ylim([0.001 1]);
xlim([0 40])