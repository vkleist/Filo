%writes all parameters into the data structure
%Parameters.<mutant>.<parametername> and saves it to the file AllParametersSimple3.mat

clear all
close all

param_set = {'WT','DLar','LiprinA','Syd1','Trio'};


for i=1:length(param_set)
    mutant = char(param_set(i));
    switch mutant
        
        case 'WT'
            % Filopodia dynamics (1/min)
            Parameters.WT.c1_sF =3.88;% short-lived filopodia birth           
            Parameters.WT.c1_LF = 1.15;%long-lived filopodia birth    
            Parameters.WT.c2_sF = 0.69;%short-lived filopodia death
            Parameters.WT.c2_LF = 0.11;%long-lived filopodia death
            %bulbous dynamics
            Parameters.WT.c3 = 0.022;% filopodia to bulbous transition 
            Parameters.WT.B50 = 0.0282;% feedback parameter 
            Parameters.WT.c4 = 1/120;% death of transient bulbous
            Parameters.WT.c5 = 0.1133;% transient to synaptogenic bulbous transition
            % time dependent functions
            % a) dampening of filo dynamics f_F(t)
            Parameters.WT.p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
            % b) increase in propensity to form bulbous f_FB(t,t_half)
            Parameters.WT.thalf = 1000;
            Parameters.WT.h = 1;
            
        case 'DLar'
            % Filopodia dynamics
            Parameters.DLar.c1_sF =2.63;%       
            Parameters.DLar.c1_LF = 1.49;%
            Parameters.DLar.c2_sF = Parameters.WT.c2_sF;%
            Parameters.DLar.c2_LF = Parameters.WT.c2_LF;%
            %bulbous dynamics
            Parameters.DLar.c3 = 0.0072;%
            Parameters.DLar.B50 = 1e-14;%
            Parameters.DLar.c4 = 0.3;%
            Parameters.DLar.c5 = 0.0228;%
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.DLar.p = Parameters.WT.p;
            % b) increase in propensity to form bulbous
            Parameters.DLar.thalf = Parameters.WT.thalf;
            Parameters.DLar.h = Parameters.WT.h;

        case 'LiprinA'   
            % Filopodia dynamics
            Parameters.LiprinA.c1_sF =3.12;%       
            Parameters.LiprinA.c1_LF = 0.99;%
            Parameters.LiprinA.c2_sF = Parameters.WT.c2_sF;%
            Parameters.LiprinA.c2_LF = Parameters.WT.c2_LF;%
            %bulbous dynamics
            Parameters.LiprinA.c3 = 0.0152;%
            Parameters.LiprinA.B50 = 0.3626;%
            Parameters.LiprinA.c4 = 0.111;%
            Parameters.LiprinA.c5 = 0.0028;%
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.LiprinA.p = Parameters.WT.p;
            % b) increase in propensity to form bulbous
            Parameters.LiprinA.thalf = Parameters.WT.thalf;
            Parameters.LiprinA.h = Parameters.WT.h;
   
        case 'Syd1'
             retract = 0;%
            % Filopodia dynamics
            Parameters.Syd1.c1_sF =2.84;%         
            Parameters.Syd1.c1_LF = 1.07;%
            Parameters.Syd1.c2_sF = Parameters.WT.c2_sF;%
            Parameters.Syd1.c2_LF = Parameters.WT.c2_LF;%
            %bulbous dynamics
            Parameters.Syd1.c3 = 0.0321;%
            Parameters.Syd1.B50 = 1.0835;%
            Parameters.Syd1.c4 = 0.169;%
            Parameters.Syd1.c5 = 0.0048;%
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.Syd1.p =  Parameters.WT.p;
            % b) increase in propensity to form bulbous
            Parameters.Syd1.thalf = Parameters.WT.thalf;
            Parameters.Syd1.h = Parameters.WT.h;
            
        case 'Trio'
             retract = 0;%0
            % Filopodia dynamics
            Parameters.Trio.c1_sF =4.71;%       
            Parameters.Trio.c1_LF = 1.61;%
            Parameters.Trio.c2_sF = Parameters.WT.c2_sF;%
            Parameters.Trio.c2_LF = Parameters.WT.c2_LF;%
            %bulbous dynamics
            Parameters.Trio.c3 = 0.0139;%
            Parameters.Trio.B50 = 0.0231;%
            Parameters.Trio.c4 = 0.1311;%
            Parameters.Trio.c5 = 0.1865;%
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.Trio.p = Parameters.WT.p;
            % b) increase in propensity to form bulbous
            Parameters.Trio.thalf = Parameters.WT.thalf;
            Parameters.Trio.h = Parameters.WT.h;
    end
end

try
    delete 'AllParametersSimple3.mat'
catch
end

save 'AllParametersSimple3.mat' Parameters 


