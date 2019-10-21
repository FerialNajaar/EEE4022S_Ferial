clear all;
close all;

%% Parameters
N = 24; 
PFA = 10^-3;   
Iterations = 500;
j = 1i;

Targets = 0;          % Number of Targets

Window_Size = 6; 
N = 32;
g_cells = 0;
g_cells1 = 0;
index = ceil((3/4)*(Window_Size*2));

SNR_dB = 20;
SNR_Linear = 10^(SNR_dB/10);

v = 10;  % Variance 

alpha_values = 1:1:500;
Pfa_values = factorial(2*Window_Size)*factorial(alpha_values+2*Window_Size-index)./(factorial(2*Window_Size-index)*factorial(alpha_values+2*Window_Size));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*PFA));
alpha = alpha_values(ind);

%% Noise Matrix 

Noise1 = zeros(N, Iterations/2);
Noise2 = zeros(N, Iterations/2);

for n = 1:N
    I_Noise = randn(1, Iterations/2);
    Q_Noise = randn(1, Iterations/2); 
%     Noise(n,:) = (I_Noise + j*Q_Noise) * (1/sqrt(2));
    Noise1(n,1:Iterations/2) = (I_Noise + j*Q_Noise) * (1/sqrt(2));
    Noise2(n, 1:Iterations/2) = sqrt(v) * (I_Noise + j*Q_Noise) * (1/sqrt(2)); % scale second half to create clutter edge
end 

Noise = [Noise1 Noise2];
%% Target 

I_Target = randn(1, Targets);
Q_Target = randn(1, Targets); 

Tgt_Voltage = sqrt(SNR_Linear); 

Tgt_Signal = Tgt_Voltage*(I_Target+ j*Q_Target)*(1/sqrt(2));

Signal = Noise;

for y = 0:(Targets -1)
    Signal(N/2,(Iterations/2) + 1 + 4*y) = Signal(N/2, (Iterations/2) + 1 + 4*y) + Tgt_Signal(y + 1);       % Single or multiple targets
%     Signal(N/2,(Iterations/2:(Iterations/2+y))) = Signal(N/2, (Iterations/2:(Iterations/2+y))) + Tgt_Signal(y + 1);   
end

Signal_Power = abs(Signal).^2; 

%% Threshold Calculation 

[num_rows, num_columns] = size(Signal_Power);     % Number of rows and columns in signal matrix 
% threshold = zeros(num_rows, num_columns);
threshold = [];

Detections_x1 = [];      % Position of Detection in row  
Detections_y1 = [];      % Position of Detection in column

% g = zeros(2*Window_Size, Iterations);
g = [];

for t = 1:num_columns                       % Iterate through columns 
    B = Signal_Power(1:num_rows,t);               % Extract a single column 

%     tempThresh = [];
    tempG = [];
    
    for i = Window_Size+(g_cells/2)+1:num_rows-Window_Size-(g_cells/2)  % Iterate through rows in extracted column
        
        CUT_Power = B(i);                   % Set the CUT 
        Window_Back = B(i-Window_Size-(g_cells/2):i-(g_cells/2)-1);     % Lagging Window
        Window_Front = B(i+1+(g_cells/2):i+(g_cells/2)+Window_Size);    % Leading Window
        
        Ref_Cells = [Window_Back; Window_Front];
        Sort_Ref_Cells = sort(Ref_Cells); 
        
        tempG = [tempG; Sort_Ref_Cells(index)];
%         g = Sum_Ref_Cells./N;          % kth value of sorted window
%         T = alpha.*g;                        % Threshold 
    
%         tempThresh = [tempThresh; T];
%         threshold(i,t) = T;
        
%         if (T<CUT_Power)
%              Detections_x1 = [Detections_x1; i];      % Add row position of detection 
%              Detections_y1 = [Detections_y1; t];      % Add column position of detection 
%         end

    end
    g = [g, tempG];
%     threshold = [threshold, tempThresh];
end

T = alpha.*g;  

%% Threshold Calculation for Ref Window sliding through rows 

[num_rows, num_columns] = size(Signal_Power);     % Number of rows and columns in signal matrix 

threshold = [];

Detections_x1 = [];      % Position of Detection in row  
Detections_y1 = [];      % Position of Detection in column

g1 = [];

for t1 = 1:num_rows                       % Iterate through columns 
    B1 = Signal_Power(t1, 1:num_columns);               % Extract a single column 

%     tempThresh = [];
    tempG1 = [];
    
    for i1 = Window_Size+(g_cells1/2)+1:num_columns-Window_Size-(g_cells1/2)  % Iterate through rows in extracted column
        
        CUT_Power = B1(i1);                   % Set the CUT 
        Window_Back1 = B1(i1-Window_Size-(g_cells1/2):i1-(g_cells1/2)-1);     % Lagging Window
        Window_Front1 = B1(i1+1+(g_cells1/2):i1+(g_cells1/2)+Window_Size);    % Leading Window
        
        Ref_Cells1 = [Window_Back1 Window_Front1];
        Sort_Ref_Cells1 = sort(Ref_Cells1); 
        
        tempG1 = [tempG1, Sort_Ref_Cells1(index)];
%         g = Sum_Ref_Cells./N;          % kth value of sorted window
%         T = alpha.*g;                        % Threshold 
    
%         tempThresh = [tempThresh; T];
%         threshold(i,t) = T;
        
%         if (T<CUT_Power)
%              Detections_x1 = [Detections_x1; i];      % Add row position of detection 
%              Detections_y1 = [Detections_y1; t];      % Add column position of detection 
%         end

    end
    g1 = [g1; tempG1];
%     threshold = [threshold, tempThresh];
end

T1 = alpha.*g1;

SP_Plot = 10*log10(Signal_Power(N/2,:));
T_Plot = 10*log10(T(N/2,:));
figure(1) 
n = 0:Iterations-1;
plot(n, SP_Plot, n, T_Plot);
ylim([-40 40]);
xlabel('Resolution Bin Index');
ylabel('Power (dB)');
title('OS-CFAR Detector: Window Sliding Through Columns');
legend('Signal','Threshold')

SP_Plot2 = 10*log10(Signal_Power(N/2,Window_Size+1:Iterations-Window_Size-g_cells1));
T_Plot2 = 10*log10(T1(N/2,:));
figure(2) 
n = 0:Iterations-(2*Window_Size+g_cells1+1);
plot(n, SP_Plot2, n, T_Plot2);
ylim([-40 40]);
xlabel('Resolution Bin Index');
ylabel('Power (dB)');
title('OS-CFAR Detector: Window Sliding Through Rows');
legend('Signal','Threshold')

