% Ferial Najaar
% NJRFER001 
% OS-CFAR Detection Algorithm

function [Detections_x, Detections_y, Signal, threshold] = OS_CFAR_Detection_2(file, PFA, g_cells,Window_Size, index, cpi, overlapFactor)
% clear all;
% file = 'AudiA1_P15_Driving_Away__20KPH_002.wav';
[dti1,speed,cpi,time] = cantenna_dop_v3_yunus(file, cpi, overlapFactor );


dti = dti1;        
% LL = dti*(1/sqrt(2));   % Linear Law Detector
LL = dti;              % Yunus: remove *1/(sqrt(2));
Signal = abs(LL).^2;    % Square Law Detector 

Target = 0; 

N = 2*(Window_Size);    % Reference Window Size 
% index = (3/4)*N;        % Ferial code 

 
Detections_x = [];      % Position of Detection in row  
Detections_y = [];      % Position of Detection in column

alpha_values = 1:1:500;
Pfa_values = factorial(N)*factorial(alpha_values+N-index)./(factorial(N-index)*factorial(alpha_values+N));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*PFA));
alpha = alpha_values(ind);

[num_rows, num_columns] = size(Signal);     % Number of rows and columns in signal matrix 
threshold = [];

for t = 1:num_columns                       % Iterate through columns 
    B = Signal(1:num_rows,t);               % Extract a single column 

    for i = Window_Size+(g_cells/2)+1:num_rows-Window_Size-(g_cells/2)  % Iterate through rows in extracted column
        CUT_Power = B(i);                   % Set the CUT 
        Window_Back = B(i-Window_Size-(g_cells/2):i-(g_cells/2)-1);     % Lagging Window
        Window_Front = B(i+1+(g_cells/2):i+(g_cells/2)+Window_Size);    % Leading Window

        Ref_Cells = [Window_Back; Window_Front];
        Sort_Ref_Cells = sort(Ref_Cells);   % Sort Reference window in ascending order 

        g = Sort_Ref_Cells(index);          % kth value of sorted window
        T = alpha*g;                        % Threshold 
    
        threshold = [threshold; T];
    
        if (T<CUT_Power)
             Detections_x = [Detections_x; i];      % Add row position of detection 
             Detections_y = [Detections_y; t];      % Add column position of detection 
             Target = Target + 1;                   % Count number of targets 
        end

    end 
end


% for t = 1:num_columns
% 
%     for i = Window_Size+(g_cells/2)+1:num_rows-Window_Size-(g_cells/2)
%         CUT_Power = Signal(i,t);
%         Window_Back = Signal(i-Window_Size-(g_cells/2):i-(g_cells/2)-1, t);
%         Window_Front = Signal(i+1+(g_cells/2):i+(g_cells/2)+Window_Size);
% 
%         Ref_Cells = [Window_Back; Window_Front];
%         Sort_Ref_Cells = sort(Ref_Cells);
% 
%         g = Sort_Ref_Cells(index);
%         T = alpha*g;
% 
%         threshold = [threshold; T];
%     
%         if (T<CUT_Power)
%              Detections_x = [Detections_x; i];
%              Detections_y = [Detections_y; t];
%              Target = Target + 1;
%         end
% 
%     end 
% end 

close all;
