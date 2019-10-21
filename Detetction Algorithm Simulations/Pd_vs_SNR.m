clear all; 
close all; 

%% Parameters

N = 24; 
PFA = 10^-6;   
Iterations = 10^6;
j = 1i;

SNR_dB_Range = 0:0.5:30;
SNR_Linear = 10.^(SNR_dB_Range/10);

Index_OS = ceil((3/4)*N);
Index_OSGO = ceil((5/12)*N);

I_Noise = randn(1, Iterations); 
Q_Noise = randn(1, Iterations); 

H0 = (I_Noise + j*Q_Noise)*(1/sqrt(2));

%% Neyman Pearson 

NP_FA = 0;
Pd_NP = [];

%% CA-CFAR 

CA_alpha = N*((PFA)^(-1/N)-1);
CA_FA = 0;
Pd_CA = [];

%% OS-CFAR
% 
% alpha_values = 1:1:500;
% Pfa_values = factorial(N)*factorial(alpha_values+N-Index_OS)./(factorial(N-Index_OS)*factorial(alpha_values+N));
% [val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*PFA));
% OS_alpha = alpha_values(ind);

OS_alpha = 16.2933;

OS_FA = 0; 
Pd_OS = [];

%% OSGO-CFAR

OSGO_alpha = 11.65;

OSGO_FA = 0; 
Pd_OSGO = [];

%% Noise 

I_Lagging = randn(N/2, Iterations);
Q_Lagging = randn(N/2, Iterations);

I_Leading = randn(N/2, Iterations); 
Q_Leading = randn(N/2, Iterations); 

Lagging_Window = (I_Lagging + j*Q_Lagging)*(1/sqrt(2));
Leading_Window = (I_Leading + j*Q_Leading)*(1/sqrt(2));

Lagging_Power = real(Lagging_Window).^2 + imag(Lagging_Window).^2;
Leading_Power = real(Leading_Window).^2 + imag(Leading_Window).^2;

Avg_Lagging = mean(Lagging_Power, 1); 
Avg_Leading = mean(Leading_Power, 1);

%% Set Thresholds

% NP Threshold 
NP_T = -log(PFA);

% CA-CFAR Threshold 
CA_g = (Avg_Lagging + Avg_Leading)/2; 
CA_T = CA_alpha.*CA_g; 

% OS-CFAR Threshold 
Ref_Cells = [Lagging_Power.', Leading_Power.'].';
OS_Sort_Ref_Cells = sort(Ref_Cells, 1);
OS_g = OS_Sort_Ref_Cells(Index_OS,:);

OS_T = OS_alpha.*OS_g;

% OSGO-CFAR Threshold
OSGO_Sort_Lag = sort(Lagging_Power,1);
OSGO_Sort_Lead = sort(Leading_Power,1);

OSGO_diff = OSGO_Sort_Lead(Index_OSGO,:) - OSGO_Sort_Lag(Index_OSGO,:);
OSGO_Greater_Lag = find(OSGO_diff < 0);
OSGO_Greater_Lead = find(OSGO_diff >= 0);

OSGO_g(OSGO_Greater_Lead) = OSGO_Sort_Lead(Index_OSGO, OSGO_Greater_Lead);
OSGO_g(OSGO_Greater_Lag) = OSGO_Sort_Lag(Index_OSGO, OSGO_Greater_Lag);

OSGO_T = OSGO_alpha.*OSGO_g;

%% Determind Pd

I_Target = randn(1,Iterations);
Q_Target = randn(1,Iterations);

for SNR = 1:length(SNR_Linear)
    
    H1 = sqrt(SNR_Linear(SNR))*(I_Target + j*Q_Target)*(1/sqrt(2)) + H0;
    H1_Power = real(H1).^2 + imag(H1).^2;
    
    % NP Detections 
    Detections = sum((H1_Power - NP_T)>0);
    Pd_NP = [Pd_NP, Detections/Iterations];
    
    % CA-CFAR Detections 
    Detections = sum((H1_Power - CA_T)>0);
    Pd_CA = [Pd_CA, Detections/Iterations];
    
    % OS-CFAR Detections
    Detections = sum((H1_Power - OS_T)>0);
    Pd_OS = [Pd_OS, Detections/Iterations];

    % OSGO-CFAR Detections
    Detections = sum((H1_Power - OSGO_T)>0);
    Pd_OSGO = [Pd_OSGO, Detections/Iterations];
    
    
end 

figure(1)
plot(SNR_dB_Range, Pd_NP, SNR_dB_Range, Pd_CA, SNR_dB_Range, Pd_OS, SNR_dB_Range, Pd_OSGO,'LineWidth',1.5);
hold on;
SNR_vs_PD_legend = legend('NP','CA-CFAR', 'OS-CFAR', 'OSGO-CFAR');
set(SNR_vs_PD_legend,'Location','SouthEast');
title('ROC for NP, CA-CFAR, OS-CFAR and OSGO-CFAR');
xlabel('SINR (dB)');
ylabel('Probability of Detection');
grid on;

figure(2)
plot(SNR_dB_Range, Pd_NP, SNR_dB_Range, Pd_CA, 'LineWidth',1.5);
hold on;
SNR_vs_PD_legend = legend('NP','CA-CFAR');
set(SNR_vs_PD_legend,'Location','SouthEast');
title('ROC for NP and CA-CFAR');
xlabel('SINR (dB)');
ylabel('Probability of Detection');
grid on;

figure(3)
plot(SNR_dB_Range, Pd_NP, SNR_dB_Range, Pd_OS, 'LineWidth',1.5);
hold on;
SNR_vs_PD_legend = legend('NP', 'OS-CFAR');
set(SNR_vs_PD_legend,'Location','SouthEast');
title('ROC for NP and OS-CFAR');
xlabel('SINR (dB)');
ylabel('Probability of Detection');
grid on;

figure(4)
plot(SNR_dB_Range, Pd_NP, SNR_dB_Range, Pd_OSGO,'LineWidth',1.5);
hold on;
SNR_vs_PD_legend = legend('NP', 'OSGO-CFAR');
set(SNR_vs_PD_legend,'Location','SouthEast');
title('ROC for NP and OSGO-CFAR');
xlabel('SINR (dB)');
ylabel('Probability of Detection');
grid on;

