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

%% alpha

PFA_error = inf;

for a_i = 0:0.05:100
    PFA_temp = 0;
    
    for m = 0:((N/2) - index)
        for p = 0:((N/2) - index)
            PFA_temp = PFA_temp + (factorial(((N/2) - index))./(factorial(m).*factorial((N/2) - index - m))).*...
                (factorial(((N/2) - index))./(factorial(p).*factorial((N/2) - index - p))).*...
                ((-1).^(N - (2*index) - m - p)./((N/2) - p)).*...
                (gamma(N - m - p).*gamma(a_i + 1)./gamma(N - m - p + a_i + 1));
        end
    end
    
    PFA_temp = (2*(index.^2)*(factorial(N/2)./(factorial((N/2) - index).*factorial(index))).^2).*PFA_temp;
    
    if abs(PFA_temp - PFA) < PFA_error
        alpha = a_i;
        PFA_error = abs(PFA_temp - PFA);
    end
end


[num_rows, num_columns] = size(Signal);     % Number of rows and columns in signal matrix 
threshold = [];

for t = 1:num_columns                       % Iterate through columns 
    B = Signal(1:num_rows,t);               % Extract a single column 

    for i = Window_Size+(g_cells/2)+1:num_rows-Window_Size-(g_cells/2)  % Iterate through rows in extracted column
        CUT_Power = B(i);                   % Set the CUT 
        Window_Back = B(i-Window_Size-(g_cells/2):i-(g_cells/2)-1);     % Lagging Window
        Window_Front = B(i+1+(g_cells/2):i+(g_cells/2)+Window_Size);    % Leading Window

        Sort_Lagging = sort(Window_Back);
        Sort_Leading = sort(Window_Front);
        
        if Sort_Leading(index) > Sort_Lagging(index)
            g = Sort_Leading(index);
           
        else
            g = Sort_Lagging(index);
        end

        T = alpha*g;                        % Threshold 
    
        threshold = [threshold; T];
    
        if (T<CUT_Power)
             Detections_x = [Detections_x; i];      % Add row position of detection 
             Detections_y = [Detections_y; t];      % Add column position of detection 
             Target = Target + 1;                   % Count number of targets 
        end

    end 
end


close all;
