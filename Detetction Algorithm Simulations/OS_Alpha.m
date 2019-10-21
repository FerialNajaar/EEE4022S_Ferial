clear all; 
close all; 

Window_Size = 12; 
N = 24; 
PFA = 10^-6;   
Iterations = 10^6;
j = 1i;
index = ceil((3/4)*N);

PFA_error = inf;

for a_temp = 0:0.0001:25
    PFA_temp = 1;
    for i = 0:(index - 1)
        PFA_temp = PFA_temp*(N - i)/(N - i + a_temp);
    end
    if abs(PFA_temp - PFA) < PFA_error
        a = a_temp;
        PFA_error = abs(PFA_temp - PFA);
    end
end