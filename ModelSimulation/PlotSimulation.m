function PlotSimulation(mutant)
%This function plots the ensemble statistics of from the Markov Model
%simulations (solid red lines are medians, dotted black lines means, light- 
% and grey areas denote the quartiles and 5-95 percentile range). 
% Input: mutant with possible entries 'WT', 'DLar', 'LiprinA', 'Trio' or 'Syd1'
% Requires that the simulations have been conducted such that the simulation files
% can be found in ./EnsembleData/

load(strcat('./EnsembleData/Simulation_',mutant))

Ref_yellow = [0.9 0.8 0.2];%

ens_filo = ens_Data(:,:,1);
ens_sFil = ens_Data(:,:,2);
ens_sbul = ens_Data(:,:,3);
ens_Lbul = ens_Data(:,:,4);
ens_syn = ens_Data(:,:,5);


%short lived Filopodia
figure(4)
hold on
XPositions = (1:time_points)./60+40;
TMP_5_50_95 = PrcTile(ens_filo,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
xlabel('time (hours)')
title(strcat('Short lived filopodia (',mutant,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 20])

%long-lived Filopodia
figure(5)
hold on
TMP_5_50_95 = PrcTile(ens_sFil,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
xlabel('time (hours)')
title(strcat('Long-lived filopodia (',mutant,')'))
ylabel('Number per Growthcone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 25])

%All filopodia
figure(10)
hold on
XPositions = (1:time_points)./60+40;
TMP_5_50_95 = PrcTile(ens_filo+ens_sFil,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
if ~strcmp(mutant,'WT')
    S = load('./EnsembleData/Simulation_WT.mat');
    Ref1 = S.ens_Data(:,:,1);
    Ref2 = S.ens_Data(:,:,2);
    MRef = mean(Ref1+Ref2);
    plot(XPositions,MRef,'Color',Ref_yellow,'LineWidth',3)
end
xlabel('time (hours)')
title(strcat('Number of filopodia (',mutant,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
ylim([0 40])


% Short-lived bulbs
figure(6)
hold on
TMP_5_50_95 = PrcTile(ens_sbul,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
plot(XPositions,mean(ens_sbul),'k:','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
if ~strcmp(mutant,'WT')
    S = load('./EnsembleData/Simulation_WT.mat');
    Ref = S.ens_Data(:,:,3);
    MRef = mean(Ref);
    plot(XPositions,MRef,'Color',Ref_yellow,'LineWidth',3)
end

xlabel('time (hours)')
title(strcat('Number of short-lived bulbous tips (',mutant,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
ylim([0 5])


% Long-lived bulbs
figure(7)
hold on
TMP_5_50_95 = PrcTile(ens_Lbul,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
plot(XPositions,mean(ens_Lbul),'k:','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
if ~strcmp(mutant,'WT')
    S = load('./EnsembleData/Simulation_WT.mat');
    Ref = S.ens_Data(:,:,4);
    MRef = mean(Ref);
    plot(XPositions,MRef,'Color',Ref_yellow,'LineWidth',3)
end

xlabel('time (hours)')
title(strcat('Number of synaptogenic bulbous tips (',mutant,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
ylim([0 5])

% All bulbs
figure(8)
hold on
TMP_5_50_95 = PrcTile(ens_Lbul+ens_sbul,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
plot(XPositions,mean((ens_sbul+ens_Lbul)),'k:','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
if ~strcmp(mutant,'WT')
    S = load('./EnsembleData/Simulation_WT.mat');
    Ref = S.ens_Data(:,:,4)+S.ens_Data(:,:,3);
    MRef = mean(Ref);
    plot(XPositions,MRef,'Color',Ref_yellow,'LineWidth',3)
end

xlabel('time (hours)')
title(strcat('Number of bulbous tips (',mutant,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
ylim([0 5])

%Synapses
figure(9)
hold on
TMP_5_50_95 = PrcTile(ens_syn,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
plot(XPositions,mean(ens_syn),'k:','LineWidth',3)
if ~strcmp(mutant,'WT')
    S = load('./EnsembleData/Simulation_WT.mat');
    Ref = S.ens_Data(:,:,5);
    MRef = mean(Ref);
    plot(XPositions,MRef,'Color',Ref_yellow,'LineWidth',3)
end

xlabel('time (hours)')
title(strcat('Number of synapses (',mutant,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
ylim([0 30])

end

function y = PrcTile(x,p,dim)
    %COMPUTES THE PERCENTILES
    [n,m] = size(x);
    
    if dim == 1
        x_sorted = sort(x,1);
        r = (p'/100)*n;   % location of percentile in array
        k = floor(r+0.5); % first index before r
        kp1 = k + 1;      % first index after r
        r = r - k;        % residual

        % set invalid indices to 1 or n
        k(k<1 | isnan(k)) = 1;
        kp1 = bsxfun( @min, kp1, n );

        % linear interpolation 
        y = (0.5+r).*x_sorted(kp1,:)+(0.5-r).*x_sorted(k,:);
        
    else
        x_sorted = sort(x,2);
        r = (p/100)*m; % location of percentile in array
        k = floor(r+0.5); % first index before r
        kp1 = k + 1;      % first index after r
        r = r - k;        % residual

        % set invalid indices to 1 or n
        k(k<1 | isnan(k)) = 1;
        kp1 = bsxfun( @min, kp1, m );

        % linear interpolation 
        y = (0.5+r).*x_sorted(:,kp1)+(0.5-r).*x_sorted(:,k);
        
    end

end

    

