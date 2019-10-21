% Ferial Najaar
% NJRFER001 
% Main code to run OS-CFAR detection algortihm 


clear all;
close all;

%% User inputs
file = 'Bremner_PoloGTI_Towards_45KPH_Lexus_Away_60KPH_003.wav';

% Yunus: moved Setup parameters for computing spectrogram
cpi = 0.1; % (s) coherent processing interval
overlapFactor = 3;  % Yunus: overlapFactor = 2 means 50% overlap 

% Yunus: moved Setup Detector 
PFA = 1e-3;             % Probability of False Alarm         
g_cells = 6;           % Guard Cells/2 for each side
Window_Size = 12;       % N/2
% g_cells = 4;           % Guard Cells/2 for each side
% Window_Size = 10; 
index = ceil((3/4)*(2*Window_Size));  % Yunus: ceil(3/4*N) 


%% Processing
[Detections_1, Detections_2, signal, threshold] = OS_CFAR_Detection_2(file, PFA, g_cells,Window_Size, index, cpi, overlapFactor);


% Setup constants and parameters
[dti,speed,cpi, overlapFactor, time] = cantenna_dop_v3_yunus(file, cpi, overlapFactor);


%for k = Detections_1          % Yunus: commented out
%    for l = Detections_2      % Yunus: commented out

        k = Detections_1;
        l = Detections_2;
        rowX = speed(k);
        colX = time(l);
        figure(1)
        text(colX, rowX,'X','FontSize',17);   % Plot X position - column then row
        
%% Yunus: implement a crude clustering algorithm  
        
        % Setup constants and parameters
[dti,speed,cpi, overlapFactor, time] = cantenna_dop_v3_yunus(file, cpi, overlapFactor);

        %  Remove all detections that occur once per time, because
        %  expecting multiple detections on the target
          IndxDup = find(diff(Detections_2) == 0); 
          IndxDupBool = ismember(Detections_2, Detections_2(IndxDup));
          TimeVector_Dupl = Detections_2(IndxDupBool);
          SpeedVector_Dupl = Detections_1(IndxDupBool);
          
        k = SpeedVector_Dupl;
        l = TimeVector_Dupl;
        rowX = speed(k);
        colX = time(l);
        figure(2)
        text(colX, rowX,'X','FontSize',10, 'Color', 'w');   % Plot X position - column then row
         k = 1;
         
         % Perform clustering: keep detection that are within one bin
         [UniqVal,TimeUniqIdx,~]=unique(TimeVector_Dupl);
         icount = 1;
         for countTimeV = 1:length(TimeUniqIdx)
            [IdxTime, ~] = find(TimeVector_Dupl == TimeVector_Dupl(TimeUniqIdx(countTimeV)));
            TimeVals =  TimeVector_Dupl(IdxTime);
            SpeedVals =  SpeedVector_Dupl(IdxTime);
            SpeedIdx = find( diff(SpeedVals) == 1) + 1;
        
            if length(SpeedIdx > 0)
            TimeVectorFinal(icount) = TimeVals(SpeedIdx(1));
            SpeedVectorFinal(icount) = SpeedVals(SpeedIdx(1));
            icount = icount + 1;
            end 
            
         end
         
         
         [dti,speed,cpi, overlapFactor, time] = cantenna_dop_v3_yunus(file, cpi, overlapFactor);

        k = SpeedVectorFinal;
        l = TimeVectorFinal;
        rowX = speed(k);
        colX = time(l);
        figure(3)
        text(colX, rowX,'X','FontSize',10, 'Color', 'w');   % Plot X position - column then row
         k = 1;
         
         
         % Remove outliers
         
         [ii,jj,kk]=unique(SpeedVectorFinal); OutliersIdx = jj(histc(kk,1:numel(ii))<2);
         SpeedVectorFinal(OutliersIdx) = [];
         TimeVectorFinal(OutliersIdx) = [];
      
        [dti,speed,cpi, overlapFactor, time] = cantenna_dop_v3_yunus(file, cpi, overlapFactor);

        k = SpeedVectorFinal;
        l = TimeVectorFinal;
        rowX = speed(k);
        colX = time(l);
        figure(4)
        text(colX, rowX,'X','FontSize',17);   % Plot X position - column then row
         k = 1;
         
         Ave_Speed = sum(rowX)/size(rowX,1);
         Ave_Speed_KPH = Ave_Speed*3.6
         
         Inst_Speed = max(rowX); 
         Inst_Speed_KPH = Inst_Speed*3.6
         % Smooth graph
%          fit = polyfit(colX, rowX', 2);
%          p = polyval(fit, rowX);
%          hold on
%          figure(4)
%          plot(colX,p, 'Color', 'm', 'LineWidth', 3);
           
     
    
%    end                     % Yunus: commented out
% end                        % Yunus: commented out

% % % make Doppler vs. time plot
% figure(3)
% imagesc(time, speed, dti);
% fit = polyfit(colX, rowX', 2);
% % imagesc(speed,time,dti);
% colormap(jet(256));
% caxis(max(dti(:)) + [-60 0]); % show 60 dB dynamic range
% colorbar;
% ylabel('Speed (m/sec)');
% xlabel('Time (sec)');
% axis xy;

