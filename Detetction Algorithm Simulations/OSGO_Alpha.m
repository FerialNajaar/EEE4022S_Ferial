clear all; 
close all; 

Window_Size = 12; 
N = 24; 
PFA = 10^-6;   
Iterations = 10^6;
j = 1i;
index = ceil((5/12)*N);

PFA_error = inf;

for a_i = 0:0.05:100
    PFA_temp = 0;
    
    for m = 0:((Window_Size) - index)
        for p = 0:((Window_Size) - index)
            PFA_temp = PFA_temp + (factorial(((Window_Size) - index))./(factorial(m).*factorial((Window_Size) - index - m))).*...
                (factorial(((Window_Size) - index))./(factorial(p).*factorial((Window_Size) - index - p))).*...
                ((-1).^(2*Window_Size - (2*index) - m - p)./((Window_Size) - p)).*...
                (gamma(2*Window_Size - m - p).*gamma(a_i + 1)./gamma(2*Window_Size - m - p + a_i + 1));
        end
    end
    
    PFA_temp = (2*(index.^2)*(factorial(Window_Size)./(factorial((Window_Size) - index).*factorial(index))).^2).*PFA_temp;
    
    if abs(PFA_temp - PFA) < PFA_error
        alpha = a_i;
        PFA_error = abs(PFA_temp - PFA);
    end
end

