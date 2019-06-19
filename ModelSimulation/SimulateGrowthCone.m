% Simulate growth cone: excecutes stochastic simulations of the data-driven
% Markov Model. 
% Input: mutant with possible entries 'WT', 'DLar', 'LiprinA', 'Trio' or 'Syd1'
% Saves the simulations into the File './EnsembleData/Simulation_<mutant>.mat'.

% Parameters of the model are loaded from the File 'AllParameters.mat'
% Run this code once for the wildtype, i.e. set mutant = 'WT' to generate
% reference simulations for the figures.

close all
clear all

%Gillespie ensemble setup
N_gil = 50; % number of stochastic simulations
T = 3600;   %final time (60 hours = 3600 minutes = 216000 seconds)
mutant = 'Syd1'; %possible entries: 'WT, 'DLar', 'LiprinA', 'Trio', 'Syd1'

%file name where calculated ensemble matrices are saved
file = strcat('EnsembleData/Simulation_',mutant); %file name where calculated ensemble matrices are saved
try
    cd 'EnsembleData'
    cd ..
catch
    mkdir('EnsembleData')
end
disp('Attention: rate Bulb -> Syn changed to simple first order reaction, line 220')

% Options for saving ensemble data
time_disc = 1;%take a snapshot every minute
time_points = ceil(T/time_disc);

%Initialize storage matrix of dimension: 
%(Nr. stoch simulations x nr. of storage time points x nr. of species)
ens_Data  = zeros(N_gil,time_points,5); %ensemble data filopodia

%Import Parameters
[c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, p, thalf, h,Initial_sFilo,Initial_LFilo] = getParameters(mutant);% 

% Stoichiometric matrix         
%   r1_sF   r1_ellF r2_sF r2_ellF   r3 r4   r5  r6    
S = [1      0       -1      0       0   0   0   0;%sF
     0      1       0       -1      0   0   0   0;%ellF
     0      0       0       0       1   -1  -1  0;%sB
     0      0       0       0       0   0   1   -1;%synB
     0      0       0       0       0   0   0   1];%S

%Start of stochastic simulations
for i=1:N_gil
    
    tic
    t = 0; %start time

    %Initialize filopodia/stab filo/bilbous/synapse
    sF = PoissRnD(Initial_sFilo); %short lived filopodia 
    LF = PoissRnD(Initial_LFilo); %long lived folopodia
    % Initial State
    X = [sF;LF;0;0;0];

    %start one growth cone simulation
    while t < T
        %compute protensities a(X,t)
        a = getPropensities(t,X,c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, thalf, h,mutant);
        a0 = sum(a);
        
        % --- Gillespie algorithm
        % 1. sample random waiting time tau from exponential distribution
        r1 = rand;
        tau = (1/a0)*log(1/r1);%
        t = t+tau; % Update time
        
        % 2. sample reaction
        r2 = rand;
        j = find(r2 <= cumsum(a)/a0,1);

        % 3. Update state according to stoichiometric vector of the choosen
        % reaction
        X = X + S(:,j);     
        % ---
        
        % ---Save simulations at the storage points into ens_Data 
        t_now = ceil(t);
        t_before = ceil(t-tau);%
        t_passed = t_before:time_disc:t_now-1;%contains all storage time points

        if t_now > T
            break;
        end

        %update ensemble matrices (storage)
        ens_Data(i,t_now,1) = X(1);
        ens_Data(i,t_now,2) = X(2);
        ens_Data(i,t_now,3) = X(3);
        ens_Data(i,t_now,4) = X(4);
        ens_Data(i,t_now,5) = X(5);


        if t_before <= 1
            continue;
        end

        %if t_passed isn't 0, fill in data for the time between
        aux = ens_Data(i,t_before,1);
        ens_Data(i,t_passed,1) = aux;
        %
        aux = ens_Data(i,t_before,2);
        ens_Data(i,t_passed,2) = aux;
        %
        aux = ens_Data(i,t_before,3);
        ens_Data(i,t_passed,3) = aux;
        %
        aux = ens_Data(i,t_before,4);
        ens_Data(i,t_passed,4) = aux;
        %
        aux = ens_Data(i,t_before,5);
        ens_Data(i,t_passed,5) = aux;
        %-----
    end

    toc

    %fill in last bits of ensemble matrices
    t_before = ceil(t-tau);
    t_passed = t_before:time_disc:T;

    aux = ens_Data(i,t_before,1);
    ens_Data(i,t_passed,1) = aux;
    %
    aux = ens_Data(i,t_before,2);
    ens_Data(i,t_passed,2) = aux;
    %
    aux = ens_Data(i,t_before,3);
    ens_Data(i,t_passed,3) = aux;
    %
    aux = ens_Data(i,t_before,4);
    ens_Data(i,t_passed,4) = aux;
    %
    aux = ens_Data(i,t_before,5);
    ens_Data(i,t_passed,5) = aux;
end


%TODO: include T, t_disc (tau),...
save(file,'ens_Data','time_points');
%plotGrowthConeSimSimple_BulbLifeTimeThreshold(time_points,ens_Data,mutant);


function a = getPropensities(t,X,c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, thalf, h,mutant)
    %loads the propensity functions
    a = zeros(1,8);
    a(1)  = c1_sF*time_dampening(t); %birth of short lived filo
    a(2)  = c1_LF*time_dampening(t); %birth of long lived filo
    a(3)  = c2_sF*X(1); %death of short lived filopodium
    a(4)  = c2_LF*X(2); %death of long lived filopodium
    a(6)  = c4*(X(3)); % short-lived bulbous death
    a(8)  = 1/180*X(4);% bulbous to synapse; alternative: 1/120*min(X(4),1); 
    
    if ~strcmp(mutant,'Trio')
        a(5)  = c3*(X(1)+X(2))*time_enhancing(t,thalf,h)*feedback(B50,X(4)); % Filopodium-to-bulbous transition
        a(7)  = X(3)*c5;%short-to-long lived bulb
    else
        a(5)  = c3*(X(1)+X(2))*time_enhancing(t,thalf,h); % Filopodium-to-bulbous transition
        a(7)  = X(3)*c5*feedback(B50,X(4));%short-to-long lived bulb
    end
end

%Function that influences the enhancement of bulbous birth: f_FB
function enhance = time_enhancing(t,thalf,h)
    scale = 2^(-h); 
    enhance = scale.*(1+tanh(3./thalf.*(t-thalf))).^(h);
end

%Function that influences the dampening of filopodia birth: f_F
function damp = time_dampening(t)
    p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
    damp = polyval(p,t);
    damp = max(1e-4,damp);
end

% Feedback function f_1
function f = feedback(B50,B)
    f = B50./(B50 + B);
end

%loads all parameters
function [c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, p, thalf, h,Initial_sFilo,Initial_LFilo] = getParameters(param)

load('AllParameters');

c1_sF = Parameters.(char(param)).c1_sF;%birth           
c1_LF = Parameters.(char(param)).c1_LF;%birth
c2_sF = Parameters.(char(param)).c2_sF;%death
c2_LF = Parameters.(char(param)).c2_LF;%death
%bulbous dynamics
c3 = Parameters.(char(param)).c3;%rF2B
B50 = Parameters.(char(param)).B50;%rB2S
c4 = Parameters.(char(param)).c4;%rB20
c5 = Parameters.(char(param)).c5;%rB2S
% time dependent functions
% a) dampening of filo dynamics; function f_F(t)
p = Parameters.(char(param)).p;
% b) increase in propensity to form bulbous; funvtion f_FB(t,thalf)
thalf = Parameters.(char(param)).thalf;
h = Parameters.(char(param)).h;

Initial_sFilo = c1_sF/c2_sF; %average number 
Initial_LFilo = c1_LF/c2_LF; %average number

end

% draws a poisson distributed random number
function r = PoissRnD(lambda)
    r = 0;
    j = 1;
    p = 0;
    while ~isempty(j)
        p = p - log(rand);
        t = (p < lambda);
        j = j(t);
        p = p(t);
        r = r+1;
    end
    r = r-1;
end