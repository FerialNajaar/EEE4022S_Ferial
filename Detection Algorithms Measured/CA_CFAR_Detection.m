% Ferial Najaar
% NJRFER001
% CA-CFAR Detection Algorithm

function [Detections_x1, Detections_y1, Signal, threshold] = CA_CFAR_Detection(file, PFA, g_cells,Window_Size, index, cpi, overlapFactor)
% clear all;
% file = 'AudiA1_P15_Driving_Away__20KPH_002.wav';
[dti1,speed,cpi,time] = cantenna_dop_v3_yunus(file, cpi, overlapFactor );

dti = dti1;        
Signal = abs(dti).^2;    % Square Law Detector 

Target = 0; 

N = 2*(Window_Size);    % Reference Window Size 
 
Detections_x1 = [];      % Position of Detection in row  
Detections_y1 = [];      % Position of Detection in column

alpha = N*((PFA)^(-1/N)-1);

[num_rows, num_columns] = size(Signal);     % Number of rows and columns in signal matrix 
threshold = [];

for t = 1:num_columns                       % Iterate through columns 
    B = Signal(1:num_rows,t);               % Extract a single column 

    for i = Window_Size+(g_cells/2)+1:num_rows-Window_Size-(g_cells/2)  % Iterate through rows in extracted column
        CUT_Power = B(i);                   % Set the CUT 
        Window_Back = B(i-Window_Size-(g_cells/2):i-(g_cells/2)-1);     % Lagging Window
        Window_Front = B(i+1+(g_cells/2):i+(g_cells/2)+Window_Size);    % Leading Window
        
        Sum_Ref_Cells = sum(Window_Back) + sum(Window_Front); 
       
        g = Sum_Ref_Cells./N;         
        T = alpha*g;                        % Threshold 
    
        threshold = [threshold; T];
    
        if (T<CUT_Power)
             Detections_x1 = [Detections_x1; i];      % Add row position of detection 
             Detections_y1 = [Detections_y1; t];      % Add column position of detection 
             Target = Target + 1;                   % Count number of targets 
        end

    end 
end

close all;
