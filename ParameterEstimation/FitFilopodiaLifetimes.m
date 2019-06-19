function[c2_sF,c2_ellF] = FitFilopodiaLifetimes()

% Analyzes whether one- or two filopodial comaprtments should be considered
% for modelling and fits parameters c2_sF and c2_ellF

load('../Data/FiloData/AllData');
%load('AllParameters');
time = {'P40';'P60'};
mutant = 'WT';

sF = [Data.(mutant).(char(time(1))).sF.LTimes;Data.(mutant).(char(time(2))).sF.LTimes];
ellF = [Data.(mutant).(char(time(1))).ellF.LTimes;Data.(mutant).(char(time(2))).ellF.LTimes];

edges = 1:59;
%%--------First part of the plot
[counts,edges] = histcounts([sF; ellF],edges);
N = sum(counts);

figure(1)
hold on
Data = counts./sum(counts);
eps = 1e-8;
VarData = StatsData(counts)+ eps;
histo = histogram('BinEdges',edges-0.5,'BinCounts',Data);
errorbar(edges(1:end-1),Data,sqrt(VarData),'LineStyle','none')

figure(2)
hold on
Data = counts./sum(counts);
eps = 1e-8;
VarData = StatsData(counts)+ eps;
histo = histogram('BinEdges',edges-0.5,'BinCounts',Data);
errorbar(edges(1:end-1),Data,sqrt(VarData),'LineStyle','none')

options = optimset('TolFun',rand,'TolX',1e-10,'Display','iter');

% 1 Compartment
edges2 = [0,edges-0.5];
[y,~] = fminsearch(@fit,rand,options,Data,edges2,VarData);
k = exp(y(1));
pred = targetfun(k,edges2);
figure(1)
plot(edges(1:end-1),pred,'b','LineWidth',2)
xlim([0 30])
ylim([-0.02 0.5])
xlabel('lifetime (in min)','FontSize',14)
ylabel('Probability','FontSize',14)
title('One Filopodia Compartment','FontSize',14)
%print(1,'-depsc2','./Figures/Fig1E_LifeTime1Compartment.eps')

% 2 compartments
r = rand*0.01;
[y,res] = fminsearch(@fit2,[r 10*r rand],options,Data,edges2,VarData);
k1 = exp(y(1));
k2 = exp(y(2));
k3 = 0.5*(1+sin(y(3)));

if k1 > k2
    c2_sF = k1;
    c2_ellF = k2;
else
    c2_sF = k2;
    c2_ellF = k1;
end
pred1 = targetfun(k1,edges2);
pred2 = targetfun(k2,edges2);
pred = k3.*pred1 + (1-k3).*pred2; 

figure(2)
h(1) = plot(edges(1:end-1),k3.*pred1,'k--','LineWidth',1);
h(2) = plot(edges(1:end-1),(1-k3).*pred2,'k-.','LineWidth',1);
h(3) = plot(edges(1:end-1),pred,'b','LineWidth',2);
xlim([0 30])
ylim([-0.02 0.5])
xlabel('lifetime (in min)','FontSize',14)
ylabel('Probability','FontSize',14)
title('Two Filopodia Compartments','FontSize',14)
legend(h,'long-lived filopodia','short-lived filopodia','sum')
%print(2,'-depsc2','./Figures/Fig1E_LifeTime2Compartments.eps')
end
%%

function res = fit(x,Data,t,VarData)
%only allow for positive entries
k = exp(x(1));
pred = targetfun(k,t);

res = sum(((pred-Data).^2)./VarData);
end

function pred = targetfun(k,t)
    pred = 1-exp(-t.*k);
    %binning ... first entry contains Prob(0 < Lifetime <= t1) 
    pred = (pred(2:end)-pred(1:end-1));
    %normalization to 1 [because we do not observe all possible times (t_i < inf)]
    pred = pred./sum(pred);
    %correction for undetected (we assume that the first entry with Lifetime 0.5 is unobserved)
    pred = pred(2:end)./(1-pred(1));
end

function res = fit2(x,Data,t,VarData)
%only allow for positive entries
k1 = exp(x(1));
k2 = exp(x(2));
k3 = 0.5*(1+sin(x(3))); %range 0 1

pred1 = targetfun(k1,t);
pred2 = targetfun(k2,t);
pred = k3.*pred1 + (1-k3).*pred2;

res = sum(((pred-Data).^2)./VarData);
end


function VarData = StatsData(Data)
%Greenwoods formula
    Total =sum(Data);
    aux = (1-cumsum(Data)./Total)./(cumsum(Data)./Total.*Total);
    VarData = (1-cumsum(Data)./Total).^2.*aux;
end


