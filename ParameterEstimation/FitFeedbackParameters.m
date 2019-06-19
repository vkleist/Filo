function[] = FitFeedbackParameters(mutantNr)

% Fits both a model without feedback (linear) and a model with
% autoinhibition feedback to the bulbous tip number distribution measured 
% during live-imaging at P60 by minimizing the Kullback-Leibler distance
% between the experimental and mode-predicted density functions


Mutants = {'WT';'DLar';'LiprinA';'Syd1';'Trio'};

FSP = 20; % finite state projection (maximum number of bulbs)

Colors = [ 0 0 1;
  0.9100 0.4100 0.1700;
  1 0 0;
  0 1 0;
  1 0 1];

c = Colors(mutantNr,:);

Mutants(mutantNr)

%Filopodia dynamics
                %WT      DLar    LiprinA     Syd1    Trio
sF          = [5.6      3.8     4.5         4.1     6.8];% average number short lived Filos @P40
LF          = [10       13      8.6         9.3     14];% average number long lived Filos @P40
%
t = 20*60;%P60
%
halfmax = 1000;
exponent = 1;
enhance = enhancing(t,halfmax,exponent);
%
p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
damp = dampening(p,t);
%expected number of sF @P60
Nbrs_sF = sF(mutantNr).*damp;%
%expected number of LF @P60
Nbrs_LF = LF(mutantNr).*damp;%
F = Nbrs_sF+Nbrs_LF;

%-------
DeathRate_sB= [1/120    0.3     0.111       0.169   0.1311];%
DeathRate_synB= [1/120  1/120  1/120      1/120  1/120]; % rate at which long-lived bulbous tips disappear (become synapses)


d_sB = DeathRate_sB(mutantNr);
d_synB = DeathRate_synB(mutantNr);


                    % n0    n1      n2      n3       n4     n5      n6
ExpDensity_synBulb = [0.025 0.65    0.21    0.11     0      0       0; %WT%[0.025 0.65    0.21    0.11     0      0       0; %WT
                     0.67   0.33    0       0        0      0       0; % DLar
                     0.76   0.24    0       0        0      0       0; %LipinA
                     0.47   0.53    0       0        0      0       0; % Syd1
                     0.1    0.9     0       0        0      0       0]; %Trio          
          
                    % n0    n1      n2      n3       n4     n5      n6
ExpDensity_sBulb = [0.89   0.11     0       0       0       0       0; %WT%[0.89   0.11     0       0       0       0       0; %WT
                    0.89   0.11     0.0056  0       0       0       0; % DLar
                    0.47   0.4      0.094   0.022   0.011   0       0; %LipinA
                    0.4    0.34     0.21    0.047   0.0033  0       0; % Syd1
                    0.57   0.25     0.11    0.025   0.029  0.021    0]; %Trio                                                                                                     
                    % n0    n1      n2      n3       n4     n5      n6
ExpDensity_AllBulb = [0.025 0.54    0.33    0.11    0       0       0;                                                
                      0.59  0.37    0.036   0.0019  0       0       0;                                            
                      0.37  0.36    0.26    0.023   0.00089 0       0;                                     
                      0.27  0.27    0.25    0.16    0.038   0.0052  0.00027;                                
                      0.055 0.46    0.3     0.14    0.025   0.0075  0.0026];  

% Compile data for the mutant
Data_synBulb = ExpDensity_synBulb(mutantNr,:);
Data_sBulb = ExpDensity_sBulb(mutantNr,:);
Data = [Data_sBulb , Data_synBulb];

% plot the data
edges = 0:7;
n = length(ExpDensity_synBulb)-1;
x_pos = 0:n;

%marginal distribution of synaptogenic bulbs synB
figure(1)
hold on
bar(edges(1:end-1),ExpDensity_synBulb(mutantNr,:)','EdgeColor','none','FaceColor',c);
hold on

%marginal distribution of transient bulbs sB
figure(2)
hold on
bar(edges(1:end-1),ExpDensity_sBulb(mutantNr,:)','EdgeColor','none','FaceColor',c);
hold on

%number distribution of all bulbs 
figure(3)
hold on
bar(edges(1:end-1),ExpDensity_AllBulb(mutantNr,:)','EdgeColor','none','FaceColor',c);
hold on

%% 1: Fit the Linear model (no feedback)
%prepare parameter fitting
StartParam = [rand rand];
options = optimset('Display','notify');
%fit parameters r3, c5 using unbounded optimization
x = fminsearch(@fitfun_linear,StartParam,options,Data,d_sB,d_synB,n);

%parameter transformations to implement bounds in otherwise unbounded optimization; back transformation
r3 = exp(x(1)); %parameter bound r3 >= 0
c5 = r3/2*(1+sin(x(2)));% parameter bound 0<= c5 < r3;

%Compute estimated probability densities 
%In the linear model, the stationary number distributions are Poisson distributions.    
%marginal distribution of transient bulbs sB
lambda1 = r3/(d_sB + c5);
p1 = poisspdf(x_pos,lambda1);

%marginal distribution of synaptogenic bulbs synB
lambda2 = c5./(d_synB).*lambda1;%
p2 = poisspdf(x_pos,lambda2);

%--- plot the fit
%marginal distribution of synaptogenic bulbs synB
figure(1)
stairs(x_pos-0.5,p2,'k:','Linewidth',4)

%marginal distribution of transient bulbs sB
figure(2)
stairs(x_pos-0.5,p1,'k:','Linewidth',4)

%Compute and plot the probability of all bulbs from the marginals
Conv_PDF = fliplr(p1'*p2);
statdistr = nan(1,n+1);
for i = 0:n
    statdistr(i+1) = sum(diag(Conv_PDF,n-i));
end

figure(3)
stairs(x_pos-0.5,statdistr,'k:','Linewidth',4)

%% 2: Data fitting feedback (product inhibition on f1 or f2)
%prepare parameter fitting
StartParam = [rand, rand,rand];
options = optimset('Display','iter');
%fit parameters r3/f_1, c5 and B50 using unbounded optimization; 
x = fminsearch(@fitfun,StartParam,options,Data,FSP,d_sB,d_synB,n,mutantNr);

if mutantNr~=5 %all mutants except for trio
    disp('Feedback on f1')
    %parameter back-transformation due to the implementation of parameter
    %bounds
    r3 = exp(x(1));%This parameter is actually r3/f_1; parameter bound [0 inf]
    
    %convert r3(t=P60)/f_1 to c3; 
    %r3(t) = f_1*f_FB(t,t_half)*(sF(t)+ellF(t))*c3 
    % <=> c3 = r3(t=P60)/[f_1*f_FB(t,t_half)*(sF(t=P60)+ellF(t=P60))]
    c3 = r3./(F.*enhance)
    if mutantNr == 4
        c5 = exp(x(2))% parameter bound [0 inf]
        B50 = 1+exp(x(3))%999/2*(1+sin(x(3)))% parameter bound [1 inf]%
    else
        c5 = r3/2*(1+sin(x(2)))% parameter bound [0 c3]
        B50 = exp(x(3))%parameter bound [0 inf]
    end
    %Make Generator
    L = makeGenerator(d_sB,d_synB,r3,B50,FSP,c5);
elseif mutantNr==5 %trio model; parameters are strictly speaking not fitted
    disp('Feedback on f2')
    %parameter back-transformation due to the implementation of parameter
    %bounds as above
    r3 = exp(x(1));
    %convert r3 @P60 to c3
    c3 = r3./(F.*enhance)*1.3
    r3 = c3*(F.*enhance);
    c5 = exp(x(2))% 
    B50 = 0.0213;% fixed
    %Make Generator
    L = makeGenerator2(d_sB,d_synB,r3,B50,FSP,c5);
end

%--- plot the fit
% Copute densities of the transient- and synaptogegic bulbs

% To compute the stationary distribution we take the generator and 
% solve for the eigenvector corresponding to eigenvalue 0
M = L';
[V,D] = eig(M);
idx = find(round(diag(D).*1e5)./1e5 == 0);
statdistr = V(:,idx)./sum(V(:,idx)); % normalization

%deconvolute to compute the marginal densities of transient- (sBulb) and
%synaptogenic (synBulb) Bulbs 
TupelM = Idx2Tupel(FSP);
P1 = zeros(n+1,1); % transient
P2 = zeros(n+1,1); % synaptogenic
for j = 0:n
    idx1 = TupelM(:,1) == j;
    P1(j+1) = sum(statdistr(idx1));
    idx2 = TupelM(:,2) == j;
    P2(j+1) = sum(statdistr(idx2));
end

%plot synaptogenic Bulbs
figure(1)
hold on
stairs(x_pos-0.5,P2,'-','Linewidth',4,'Color',[0.2 0.2 0.2])
legend('data','no feedback', 'feedback')

%plot transient Bulbs
figure(2)
hold on
stairs(x_pos-0.5,P1,'-','Linewidth',4,'Color',[0.2 0.2 0.2])
legend('data','no feedback', 'feedback')

%plot density of all Bulbs
Prob_allB = nan(n+1,1);
for j = 0:n
    idx = sum(TupelM,2)==j;
    Prob_allB(j+1) = sum(statdistr(idx));
end
figure(3)
stairs(x_pos-0.5,Prob_allB,'-','Linewidth',4,'Color',[0.2 0.2 0.2])
ylabel('probability','FontSize',14)
xlabel('number bulbous tips/(GC and time)','FontSize',14)
title(strcat(Mutants(mutantNr),': all bulbs'),'FontSize',16)
legend('data','no feedback', 'feedback')


%% 3: Compute Parameters of Table S1

disp('Prob. Bulb becomes synaptogene (feedback model, feedback off):')
PsB2synB = c5/(c5+d_sB)
disp('Average bulbs @P60 (predicted)')
x_pos*Prob_allB

Expect_f = 0;
for j = 0:n
    Expect_f = Expect_f + P2(j+1) * feedback(B50,j);
end
disp('average r_5 @P60')
if mutantNr == 5
    (x_pos*P1)*c5.*Expect_f
    Expect_f = 1;
else
    (x_pos*P1)*c5
end

disp('average feedback f1 @P60')
Expect_f
disp('average r_2B @P60')
Expect_r2B = r3
disp('average r_3 @P60')
r3*Expect_f
disp('average r_4 @P60')
(x_pos*P1)*d_sB
disp('Average bulbs @P60 (measured)')
x_pos*ExpDensity_AllBulb(mutantNr,:)'

end
% -----Accesory functions-----

function res = fitfun_linear(x,Data,d_sB,d_synB,n)
    %parameter transformations to implement parameter bounds
    r3 = exp(x(1)); %bound [0 inf]
    c5 = r3/2*(1+sin(x(2)));%bound [0 r3]
    
    %In the linear model, the stationary number distributions are Poisson
    %distributions.
    
    %marginal distribution of transient bulbs sB
    x_pos = 0:n;
    lambda1 = r3/(d_sB + c5);
    p1 = poisspdf(x_pos,lambda1);
    
    %marginal distribution of synaptogenic bulbs synB
    lambda2 = sum(c5.*p1./(d_synB));
    p2 = poisspdf(x_pos,lambda2);

    statdistr = [p1,p2];
    %compute Kullback-Leibler distance between data and predicted marginal
    %distributions
    res = KL(Data,statdistr);
end


function res = fitfun(x,Data,FSP,d_sB,d_synB,n,mutant)
    if mutant ~=5 %all mutants except Trio
        %Implement parameter bounds by transformation
        r3 = exp(x(1));% parameter bound [0 inf]
        if mutant == 4
            c5 = exp(x(2));% parameter bound [0 inf]
            B50 = 1+exp(x(3));% parameter bound [1 inf]%
        else
            c5 = r3/2*(1+sin(x(2)));% parameter bound [0 c3]
            B50 = exp(x(3));% parameter bound [0 inf]
        end
        L = makeGenerator(d_sB,d_synB,r3,B50,FSP,c5);
    elseif mutant==5 %trio model; parameters are strictly speaking not fitted
        r3 = exp(x(1));% parameter bound [0 inf]
        c5 = exp(x(2));% parameter bound [0 inf]
        B50 = 1e-4;%
        %Make Generator
        L = makeGenerator2(d_sB,d_synB,r3,B50,FSP,c5);
    end
    
    %solve for stationary distribution
    M = L';
    [V,D] = eig(M);
    %find eigenvecor corresponding to eigenvalue 0 
    idx = find(round(diag(D).*1e5)./1e5 == 0);
    statdistr = V(:,idx)./sum(V(:,idx));
    
    % Extract marginal distributions for sBulb (tansient) and synBulb (synaptogenic) from the
    % stationary distribution
    
    %converts the indices in the generator to tupels
    TupelM = Idx2Tupel(FSP);
    
    %Compute the number distribution for transient- and synaptogenic bulbs
    P1 = zeros(1,n+1);
    P2 = zeros(1,n+1);
    for j = 0:n
        idx1 = TupelM(:,1) == j;
        P1(j+1) = sum(statdistr(idx1));%marginal distribution of s_Bulb
        idx2 = TupelM(:,2) == j;
        P2(j+1) = sum(statdistr(idx2));%marginal distribution of syn_Bulb
    end
    %compute Kullback-Leibler distance between data and predicted marginal
    %distributions
    res = KL(Data,[P1, P2]);
end

function L = makeGenerator(d_sB,d_synB,r3,B50,FSP,c5)
    %makes the infinitesimal generator for all mutants except trio

    L = zeros((FSP+1)^2,(FSP+1)^2);
    
    for i = 0:FSP % number short lived bulbs
        for j = 0:FSP %number long lived bulbs
            idx = i*(FSP+1)+(j+1);%current state
            
            if i > 0 %reaction: sB -> *
                idx_death_i = (i-1)*(FSP+1)+(j+1);%(i-1,j)
                L(idx,idx_death_i) = d_sB*i;
            end
            
            if j > 0 %reaction:synB -> *
                idx_death_j = (i)*(FSP+1)+(j);%(i,j-1)
                L(idx,idx_death_j) = (d_synB)*j;
            end
            
            if i < FSP %reaction: F -> sB 
                idx_birth_i = (i+1)*(FSP+1)+(j+1);%(i+1,j)
                L(idx,idx_birth_i) = r3.*feedback(B50,j);%
            end
            
            if j < FSP && i > 0%reaction: sB -> synB
                idx_birth_j = (i-1)*(FSP+1)+(j+2);%(i-1,j+1)
                L(idx,idx_birth_j) = (i)*c5;%
            end

        end
    end
    L = L - diag(sum(L,2)); %make it a generator matrix (row sum 0)
end

function L = makeGenerator2(d_sB,d_synB,r3,B50,FSP,c5)
    %makes the infinitesimal generator for all mutants except trio
    L = zeros((FSP+1)^2,(FSP+1)^2);
    
    for i = 0:FSP % number short lived bulb
        for j = 0:FSP %number long lived buld
            idx = i*(FSP+1)+(j+1);%current state
            
            if i > 0 %reaction:sB -> *
                idx_death_i = (i-1)*(FSP+1)+(j+1);%(i-1,j)
                L(idx,idx_death_i) = d_sB*i;
            end
            
            if j > 0 %reaction:synB -> *
                idx_death_j = (i)*(FSP+1)+(j);%(i,j-1)
                L(idx,idx_death_j) = (d_synB)*j;
            end
            
            if i < FSP %reaction:F -> sB +1
                idx_birth_i = (i+1)*(FSP+1)+(j+1);%(i+1,j)
                L(idx,idx_birth_i) = r3;%.*feedback(B50,(i+j));
            end
            
            if j < FSP && i > 0%reaction: sB -> synB
                idx_birth_j = (i-1)*(FSP+1)+(j+2);%(i-1,j+1)
                L(idx,idx_birth_j) = (i)*c5*feedback(B50,(j));%
            end

        end
    end
    L = L - diag(sum(L,2));%make it a generator matrix (row sum 0)
end


function f = feedback(B50,B)
%feedback function for the auto-inhibition (product inhibition)
    f = B50./(B50 + B);
end

function p = poisspdf(x,lambda)
    %Poisson probability density function (pdf)
    p = lambda.^x./factorial(x).*exp(-lambda); 
end

function TupelM = Idx2Tupel(n)
    TupelM = nan((n+1)^2,2);
    %This Matrix contains at row i the number of short-lived bulbs (first
    %column) and long lived bulbs (second column)
    counter = 0;
    for i = 0:n
        for j = 0:n
            counter = counter +1;
            TupelM(counter,1) = i;
            TupelM(counter,2) = j;
        end
    end
end

function enhance = enhancing(t,thalf,h)
    % function f_FB(t,t_half)
    scale = 2^(-h); 
    enhance = scale.*(1+tanh(3./thalf.*(t-thalf))).^(h);
end

function d = dampening(p,t)
    %function f_F(t)
    d = polyval(p,t);
    d = max(1e-4,d);
end

function D = KL(p1,p2)
    %Computes the kullback-Leibler distance between two probability densities
    %on a discrete state space.
    % p1 and p2 are probabilities in some discrete state space (pdf)
    p1(p1==0) = 1e-13;
    p2(p2==0) = 1e-13;

    D = sum(p1.*log(p1./p2));
end

