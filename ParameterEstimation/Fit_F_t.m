function [] = Fit_F_t()
% Fits the function f_F(t) (a polynomial) to filopodia numbers from fixed preparations at
% P40, P50, ..., P100


format short e
% clear all
% close all

eps = 1e-6;

%Data [number of filopodia] from fixed preparations
Data   = [20.9 13.6   12.3    10.8    3.6       0.4       0];
StdData = ([4.7  3.6    3.2     2.8     1.7       0.8       eps]);
%normalize data and variance
CVData = StdData./Data;
Data = Data./Data(1); 
StdData = CVData.*Data;

t =    [0 10*60 20*60 30*60 40*60   50*60   60*60];

%1. Fit polynomial
figure(2)
hold on
plot(t./60+40, Data,'ko','MarkerSize',14,'MarkerFaceColor','k')
tplot = linspace(0,t(end));
for n = 1:5%length(Data)
    [p,S] = polyfit(t,Data,n);
    pred = polyval(p,tplot);
    h(n) = plot(tplot./60+40, pred);
end
xlim([39 101])
title('polynome')
xlabel('time [hours]')
ylabel('normalized filopodia count')
legend(h,num2str((1:5)'))
p

%plot time-dependent functions
figure(1)
hold on
errorbar(t./60+40,Data,StdData)
plot(t./60+40, Data,'ko','MarkerSize',14,'MarkerFaceColor','k')
damp = pred;
damp = max(1e-6,damp);

%the function f_FB(t,t_half) for the enhancement of the
%filopodia-to-bulbous transition
halfmax = 1000;
exponent = 1;
scale = 2^(-exponent);
damp = polyval(p,tplot);
damp = max(1e-6,damp);
figure(1)
h(1) = plot(tplot./60+40, damp,'r--','Linewidth',3);
enhance = scale*(1+tanh(3/halfmax*(tplot-halfmax))).^(exponent);
h(2) = plot(tplot./60+40, enhance,'k--','Linewidth',4);
xlabel('time (hours)','Fontsize',16)
legend(h,'Filopodial birth rate','Rate of bulbous development')
ylabel('Reaction rate (relative to maximum)','Fontsize',16)
title('Time-dependent reactions','Fontsize',18)
set(gca,'FontSize',16);



