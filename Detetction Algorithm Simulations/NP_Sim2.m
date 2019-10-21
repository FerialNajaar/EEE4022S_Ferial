clear all;
close all;

%% Constant Parameters 
N = 1;  
PFA = 10^-3;   
Iterations = 500;
j = 1i;

Targets = 1;          % Number of Targets

SNR_dB = 20;
SNR_Linear = 10^(SNR_dB/10);

p = 0:Iterations-1;

%% Noise Matrix 

Noise = zeros(N, Iterations); 

for n = 1:N
    I_Noise = randn(1, Iterations);
    Q_Noise = randn(1, Iterations); 
    Noise(n,:) = (I_Noise + j*Q_Noise) * (1/sqrt(2));
end 

SquareLawDetect = abs(Noise).^2;

%% Target 

I_Target = randn(1, Targets);
Q_Target = randn(1, Targets); 

Tgt_Voltage = sqrt(SNR_Linear); 

Tgt_Signal = Tgt_Voltage*(I_Target+ j*Q_Target)*(1/sqrt(2));

Signal = Noise;

Signal(Iterations/2) = Signal(Iterations/2) + Tgt_Signal; 

Signal_Power = abs(Signal).^2; 

%% Threshold Calculation 
T = -log(PFA);      % Threshold 
Const_T = T*ones(1, length(p));

figure(1)
n = 0:Iterations-1;
plot(n, 20*log10(SquareLawDetect));
ylim([-30 30]);
xlabel('Resolution Bin Index');
ylabel('Power (dB)');
title('Gaussian Noise Signal');

figure(2) 
n = 0:Iterations-1;
plot(n, 20*log10(Signal_Power));
ylim([-40 40]);
xlabel('Resolution Bin Index');
ylabel('Power (dB)');
title('Target in Noise Signal');

figure(3) 
n = 0:Iterations-1;
plot(n, 20*log10(Signal_Power), n, 20*log10(Const_T));
ylim([-40 40]);
xlabel('Resolution Bin Index');
ylabel('Power (dB)');
title('Neyman-Pearson Detector');