close all
clear all

% Depicts the number of filopodia 
%

Mutants = {'WT';'DLar'; 'LiprinA'; 'Syd1';'Trio'};
Times = {'P40';'P60'};

load('FiloData/AllData.mat');

MeanNbrsFilos = zeros(length(Mutants),2);
SDNbrsFilos = zeros(length(Mutants),2);
MeanNbrsStabFilos = zeros(length(Mutants),2);
SDNbrsStabFilos = zeros(length(Mutants),2);

Var1 = cell(length(Mutants),1);
Var2 = cell(length(Mutants),1);

Var3 = cell(length(Mutants),1);
Var4 = cell(length(Mutants),1);


for j = 1:length(Mutants)
    for i = 1:length(Times)    
        Exp = unique(Data.(char(Mutants(j))).(char(Times(i))).sF.GC);
        NrExp = length(Exp);
        
        NbrsStabFilo = zeros(NrExp,60);
        NbrsFilos = zeros(NrExp,60);
        
        for counter1 = 1:NrExp
            %For each growth cone
            %Number of Filopodia at time instance
            Idx = Data.(char(Mutants(j))).(char(Times(i))).sF.GC == Exp(counter1);
            NrFilos = sum(Idx);
            Starts = Data.(char(Mutants(j))).(char(Times(i))).sF.StartTimes(Idx);
            Ends = Data.(char(Mutants(j))).(char(Times(i))).sF.EndTimes(Idx);
            AllNrFilos = zeros(1,60);
            for counter2 = 1:NrFilos
                StartT = Starts(counter2)+1;
                EndT = Ends(counter2)+1;
                AllNrFilos(StartT:EndT) = AllNrFilos(StartT:EndT)+1;
            end
            % %Number of Filopodia at time instance
            Idx = Data.(char(Mutants(j))).(char(Times(i))).ellF.GC == Exp(counter1);
            NrStabFilos = sum(Idx);
            Starts = Data.(char(Mutants(j))).(char(Times(i))).ellF.StartTimes(Idx);
            Ends = Data.(char(Mutants(j))).(char(Times(i))).ellF.EndTimes(Idx);
            AllStabFilos = zeros(1,60);   
            for counter2 = 1:NrStabFilos
                StartT = Starts(counter2)+1;
                EndT = Ends(counter2)+1;
                AllStabFilos(StartT:EndT)= AllStabFilos(StartT:EndT)+ 1;
            end
            
            %plot Number of Filopodia at time instance
            figure(101+j)
            subplot(2,3,(i-1)*3+1)
            hold on
            plot(AllNrFilos,':','LineWidth',2)
            subplot(2,3,(i-1)*3+2)
            hold on
            plot(AllStabFilos,':','LineWidth',2)
            
            NbrsFilos(counter1,:) = NbrsFilos(counter1,:) + AllNrFilos;
            NbrsStabFilo(counter1,:) = NbrsStabFilo(counter1,:) + AllStabFilos;
        end % end Exp
        
        %short lived filopodia
        data = NbrsFilos(:);
        m = mean(data);
        s = std(data);
        %save data 
        if i == 1 % P40
            Var1{j} = strcat(num2str(m,2),'(',num2str(s,2),')');
        else %P60
            Var2{j} = strcat(num2str(m,2),'(',num2str(s,2),')');
        end  
              
        
        if i == 1 && j == 1 %WT P40
            NbrsWtFilosP40 = NbrsFilos;
            NbrsWtStabFilosP40 = NbrsStabFilo;
        elseif i == 2 && j == 1 %WT P60
            NbrsWtFilosP60 = NbrsFilos;
            NbrsWtStabFilosP60 = NbrsStabFilo;
        end
        
        %Plot number histogram
        edges = (0:20)-0.5;
        figure(201+j)
        subplot(2,2,(i-1)*2+1)
        hold on
        title(strcat('Nr. (detect) filo.(',char(Mutants{j}),')'),'FontSize',16)
        [counts,edges] = histcounts(data,edges);
        histogram('BinEdges',edges,'BinCounts',counts./sum(counts))
        line([m  m],[0 0.2],'Color','r','LineStyle','--','LineWidth',3)
        line([m+s m+s],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
        line([m-s  m-s],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
        ylim([0 0.3])
        
        
        %long lived filopodia
        data = NbrsStabFilo(:);
        m = mean(data);
        s = std(data);
        %save data 
        if i == 1 % P40
            Var3{j} = strcat(num2str(m,2),'(',num2str(s,2),')');
        else %P60
            Var4{j} = strcat(num2str(m,2),'(',num2str(s,2),')');
        end  
        
        subplot(2,2,(i-1)*2+2)
        hold on
        title(strcat('Nr. stab. filo.(',char(Mutants{j}),')'),'FontSize',16)
        [counts,edges] = histcounts(data,edges);
        histogram('BinEdges',edges,'BinCounts',counts./sum(counts))
        line([m  m],[0 0.2],'Color','r','LineStyle','--','LineWidth',3)
        line([m+s m+s],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
        line([m-s  m-s],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
        ylim([0 0.3])
        
        
        %Plot mean number of Filopodia at time instance
        figure(101+j)
        subplot(2,3,(i-1)*3+1)
        hold on
        plot(mean(NbrsFilos),'k-','LineWidth',3)
        subplot(2,3,(i-1)*3+2)
        hold on
        plot(mean(NbrsStabFilo),'k-','LineWidth',3)
        
        %Plot ratios of ellF:sF at time instance
        subplot(2,3,(i-1)*3+3)
        hold on
        ratio = mean(NbrsFilos)./mean(NbrsStabFilo);
        plot(ratio,'k-','LineWidth',3)  
        ratio(isnan(ratio)) = inf;
        ratio = ratio(~isinf(ratio));
        line([1 60],[mean(ratio) mean(ratio)],'LineStyle',':','Color','r','LineWidth',3)
        ylim([0 3])
        title('Ratio: Filo-to-stab.Filo')

        
        Data.(char(Mutants(j))).(char(Times(i))).sF.Numbers = NbrsFilos(:);
        Data.(char(Mutants(j))).(char(Times(i))).ellF.Numbers = NbrsStabFilo(:);
    end % end over times i
    
    figure(101+j)
    subplot(2,3,1)
    ylabel('P40','FontWeight','bold','FontSize',14)
    title(strcat('Nr. (detectable) filo. (',char(Mutants{j}),')'))
    subplot(2,3,2)
    title(strcat('Nr. stab. filo.(',char(Mutants{j}),')'))
    subplot(2,3,4)
    ylabel('P60','FontWeight','bold','FontSize',14)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    subplot(2,3,5)
    xlabel('time  (min)','FontWeight','bold','FontSize',12)
    subplot(2,3,6)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    %print(101+j,'-depsc2',strcat('./Figures/FiloNumbers',char(Mutants{j}),'.eps'))
    
    figure(201+j)
    subplot(2,2,1)
    ylabel('P40','FontWeight','bold','FontSize',14)
    %title(strcat('Nr. (detectable) filo. (',char(Mutants{j}),')'))
    subplot(2,2,3)
    ylabel('P60','FontWeight','bold','FontSize',14)
    %print('-depsc2',strcat('./Figures/FiloNumberDistribution',char(Mutants{j}),'.eps'));
end % end over Mutants

%Print mean numbers into table
Names = {'Mutant';'P40';'P60'};
%print data to table
T = table(Mutants,Var1,Var2);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Average (standard deviation) numbers of short-lived filopodia per time instance';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);    


%Print mean numbers into table
Names = {'Mutant';'P40';'P60'};
%print data to table
T = table(Mutants,Var3,Var4);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Average (standard deviation) numbers of long-lived filopodia per time instance';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);    

%--------
