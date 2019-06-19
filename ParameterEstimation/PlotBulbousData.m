function [] = PlotBulbousData()

%plots the bulbous tip data (number per time instance) and generates LaTeX
%tables of the data.

LifeTimeThreshold = 40; %min if Lifetime >= threshold the bulbous tip is considered synaptogenic
Mutants = {'WT';'DLar';'LiprinA';'Syd1';'Trio'};

%Lifetime of filopodia that, at some point had a bulbous tip (min); 
Data.BulbFilo.WT = [60 42   60  14 60  15  60  12  60  60  49  60  60];%lifetime
Data.GC.WT =       [1   1   2   3  3   3   4   4   5   5   6   6   6];%growth cone
Data.BulbFilo.DLar = [9 20  60  24  11 60   10];
Data.GC.DLar =         [1 1   2   3   3   3   3];
Data.BulbFilo.LiprinA = [60 12 47 17 23 23 53 60 19 13 31 17 13 43 14];
Data.GC.LiprinA =       [1  1  1  1  1  1  2  2  2  2  3  3  3  3   3];
Data.BulbFilo.Syd1 = [7 12 5 19 60 20 6 12 23 6 35 12 17  7  2  26 13 36 21 32 28 23 60 1 13 3 2 7 60 46 47 6 7 21 19 7 6 2 17 6 5 16 60 26 60 38 2 2 10];
Data.GC.Syd1 =       [1 1  1  1  1  1 1 1  1  1 1  1  1   2  2  2  3  3  3  3  3   3 3  3 3  4 4 4  4 4  4  4 4 4  4  4 4 4 4  4 4 5  5   5 5  5  5 5 5];
Data.BulbFilo.Trio = [7 60 17 15    11 21 5 13  10 13 20 60 60 31 7 4 60 60 1 3 20 11 21 18 18 60 60 57];
Data.GC.Trio =       [1 1   1 1     1   1 1 1   1  2  2  2  2  2  2 3 3  3  3 4 4   4 4  4  4   4  4 4];

%Lifetime of only the bulbous tip on the bulbous filopodia (min)
Data.BulbTip.WT = [60 28 53 4 58 8 60 1 60 60 40 60 56];
Data.BulbTip.DLar = [3 1 9 4 2 60 1];
Data.BulbTip.LiprinA = [6 2 8 9 1 6 43 32 1 2 6 8 7 29 9];
Data.BulbTip.Syd1 = [2 7 3 4 41 1 1 9 10 3 12 9 1 1 1 2 12 23 3 1 6 17 2 1 3 1 1 2 59 10 27 2 5 20 15 2 1 1 10 1 1 10 60 13 4 8 2 1 2];
Data.BulbTip.Trio = [4 59 10 4 5 15 2 8 4 1 3 60 12 18 6 2 6 53 1 2 11 3 15 3 2 25 43 21];

%Emergence of the bulbous filopodium in min after the start of live imaging
Data.BulbFiloStart.WT = [0 0 0 23 0 38 0 25 0 0 6 0 0];
Data.BulbFiloStart.DLar = [2 4 0 0 14 0 37];
Data.BulbFiloStart.LiprinA = [0 0 13 0 30 0 7 0 25 32 12 19 0 17 0];
Data.BulbFiloStart.Syd1 = [45 25 25 0 0 39 0 48 0 35 20 48 5 35 37 0 47 21 39 23 14 37 0 11 4 37 0 17 0 14 0 36 12 39 41 29 5 48 31 0 0 40 0 34 0 14 58 7 6];
Data.BulbFiloStart.Trio = [42 0 43 6 3 0 13 16 50 47 7 0 0 29 19 0 0 0 38 8 5 10 0 0 2 0 0 3];

%preparation of tables 
Var1 = cell(length(Mutants),1); % mean (std) nr bulbous tips

Var2 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var3 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var4 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var5 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var6 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var7 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var8 = cell(length(Mutants),1); % prob that there are 0,1,...5


Var2_1 = cell(length(Mutants),1); % mean (std) nr bulbous tips
Var2_2 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_3 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_4 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_5 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_6 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_7 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_8 = cell(length(Mutants),1); % prob that there are 0,1,...5

%Compute Nr BulbousTips per Time Instance
ReSample = 0;%100;%for each GrowthCone; if 0 then no resample

for j = 1:length(Mutants)
    mutant = char(Mutants(j));
    GCs = unique((Data.GC.(mutant)));
    NrGCs = length(GCs);
    
    BulbTips = Data.BulbTip.(char(Mutants(j)));
    %--------------
    %if j > 2 % mutant
        %lgx = isoutlier(BulbTips);
        lgx = BulbTips>=LifeTimeThreshold;
        if mean(BulbTips(lgx)) > mean(BulbTips(~lgx))
            Long = lgx;
            Short = ~lgx;
        else
            Long = ~lgx;
            Short = lgx;
        end

    %--------------
    
    if ReSample == 0
        Nbrs_sBulbous = zeros(NrGCs*1,60);
        Nbrs_synBulbous = zeros(NrGCs*1,60);
    else
        Nbrs_sBulbous = zeros(NrGCs*ReSample,60);
        Nbrs_synBulbous = zeros(NrGCs*ReSample,60);
    end
    ReSampleCounter = 0;
    for counter1 = GCs
        %long lived
        Idx = Data.GC.(char(Mutants(j))) == counter1 & Long;
        NrSynBulbous = sum(Idx);
        Tmp = Data.BulbFiloStart.(char(Mutants(j)))+1;
        StartsSynBF =Tmp(Idx); %Bulbous Filo
        Tmp = Data.BulbFilo.(char(Mutants(j)));
        LT_SynBF = Tmp(Idx);%LifeTime of Bulbous Filo
        Tmp = Data.BulbTip.(char(Mutants(j)));
        LT_SynBT = Tmp(Idx);%LifeTime of Bulbous Tip
        %short-lived
        Idx = Data.GC.(char(Mutants(j))) == counter1 & Short;
        

        NrSBulbous = sum(Idx);
        Tmp = Data.BulbFiloStart.(char(Mutants(j)))+1;
        StartsSBF =Tmp(Idx); %Bulbous Filo
        Tmp = Data.BulbFilo.(char(Mutants(j)));
        LT_SBF = Tmp(Idx);%LifeTime of Bulbous Filo
        Tmp = Data.BulbTip.(char(Mutants(j)));
        LT_SBT = Tmp(Idx);%LifeTime of Bulbous Tip
        
        MaxStartShift_syn = LT_SynBF-LT_SynBT;
        MaxStartShift_s = LT_SBF-LT_SBT;
        
        Starts_syn = nan(1,NrSynBulbous);
        Starts_s = nan(1,NrSBulbous);
         
        if ReSample == 0
            for counter2 = 1:NrSynBulbous
               Starts_syn(counter2) = StartsSynBF(counter2);
            end
            for counter2 = 1:NrSBulbous
               Starts_s(counter2) = StartsSBF(counter2);
            end
            Ends_syn = LT_SynBT+Starts_syn-1;%Bulbous Tip
            Ends_s = LT_SBT+Starts_s-1;%Bulbous Tip
            
            AllSynBulbous = zeros(1,60);
            for counter2 = 1:NrSynBulbous
                StartT = Starts_syn(counter2);
                EndT = Ends_syn(counter2);
                AllSynBulbous(StartT:EndT) = AllSynBulbous(StartT:EndT)+1;
            end
            
            AllSBulbous = zeros(1,60);
            for counter2 = 1:NrSBulbous
                StartT = Starts_s(counter2);
                EndT = Ends_s(counter2);
                AllSBulbous(StartT:EndT) = AllSBulbous(StartT:EndT)+1;
            end
            %plot Number of Filopodia at time instance
            figure(100+j)
            hold on
            plot(AllSynBulbous,':','LineWidth',2)
            
            figure(300+j)
            hold on
            plot(AllSBulbous,':','LineWidth',2)

            ReSampleCounter = ReSampleCounter+1;
            Nbrs_synBulbous(ReSampleCounter,:) = Nbrs_synBulbous(ReSampleCounter,:) + AllSynBulbous;
            Nbrs_sBulbous(ReSampleCounter,:) = Nbrs_sBulbous(ReSampleCounter,:) + AllSBulbous;

        else
            for z = 1:ReSample
                %Syn bulbous
                for counter2 = 1:NrSynBulbous
                    if MaxStartShift_syn(counter2) == 0
                        Starts_syn(counter2) = StartsSynBF(counter2);
                    else
                        Starts_syn(counter2) = StartsSynBF(counter2) + randi(MaxStartShift_syn(counter2)+1)-1;
                    end
                end
                
                Ends_syn = LT_SynBT+Starts_syn-1;%Bulbous Tip
                AllSynBulbous = zeros(1,60);
                for counter2 = 1:NrSynBulbous
                    StartT = Starts_syn(counter2);
                    EndT = Ends_syn(counter2);
                    AllSynBulbous(StartT:EndT) = AllSynBulbous(StartT:EndT)+1;
                end
                ReSampleCounter = ReSampleCounter+1;
                Nbrs_synBulbous(ReSampleCounter,:) = Nbrs_synBulbous(ReSampleCounter,:) + AllSynBulbous;
                
                %plot Number of Filopodia at time instance
                figure(100+j)
                hold on
                plot(AllSynBulbous,':','LineWidth',2)
                
                %Short-lived Bulbous
                for counter2 = 1:NrSBulbous
                    if MaxStartShift_s(counter2) == 0
                        Starts_s(counter2) = StartsSBF(counter2);
                    else
                        Starts_s(counter2) = StartsSBF(counter2) + randi(MaxStartShift_s(counter2)+1)-1;
                    end
                end
                
                Ends_s = LT_SBT+Starts_s-1;%Bulbous Tip
                AllSBulbous = zeros(1,60);
                for counter2 = 1:NrSBulbous
                    StartT = Starts_s(counter2);
                    EndT = Ends_s(counter2);
                    AllSBulbous(StartT:EndT) = AllSBulbous(StartT:EndT)+1;
                end
                
                                %plot Number of Filopodia at time instance
                figure(300+j)
                hold on
                plot(AllSBulbous,':','LineWidth',2)
                
                Nbrs_sBulbous(ReSampleCounter,:) = Nbrs_sBulbous(ReSampleCounter,:) + AllSBulbous;

            end
        end
    end % end over growth cone
    %short lived filopodia
    data_syn = Nbrs_synBulbous(:);
    m_syn = mean(data_syn);
    s_syn = std(data_syn);
    
   %save data 
    Var1{j} = strcat(num2str(m_syn,2),'(',num2str(s_syn,2),')'); 
    
    data_s = Nbrs_sBulbous(:);
    m_s = mean(data_s);
    s_s = std(data_s);
    %save data 
    Var2_1{j} = strcat(num2str(m_syn,2),'(',num2str(s_syn,2),')'); 

    %Plot mean number of bulbous at time instance
    figure(100+j)
    hold on
    plot(mean(Nbrs_synBulbous),'k-','LineWidth',3)
    title(strcat('Nr. synaptogenic bulbous tips (',char(Mutants{j}),')'),'FontSize',14)
    ylabel('P60','FontWeight','bold','FontSize',14)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    %print(100+j,'-depsc2',strcat('./Figures/SynBulbousNumbers',char(Mutants{j}),'.eps'))
    
    figure(300+j)
    hold on
    plot(mean(Nbrs_sBulbous),'k-','LineWidth',3)
    title(strcat('Nr. short-lived bulbous tips (',char(Mutants{j}),')'),'FontSize',14)
    ylabel('P60','FontWeight','bold','FontSize',14)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    %print(300+j,'-depsc2',strcat('./Figures/ShortBulbousNumbers',char(Mutants{j}),'.eps'))
    
    %Plot number histogram
    edges = (0:10)-0.5;
    figure(200+j)
    hold on
    title(strcat('Nr. synaptogenic bulbous  tips (',char(Mutants{j}),')'),'FontSize',14)
    [counts,edges] = histcounts(data_syn,edges);
    X = counts./sum(counts);
    histogram('BinEdges',edges,'BinCounts',X)
    line([m_syn  m_syn],[0 0.2],'Color','r','LineStyle','--','LineWidth',3)
    line([m_syn+s_syn m_syn+s_syn],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    line([m_syn-s_syn  m_syn-s_syn],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    ylim([0 0.9])
    ylabel('Frequency/Probability','FontWeight','bold','FontSize',14)
    xlabel('Number bulbous tips')
    %print(200+j,'-depsc2',strcat('./Figures/SynBulbousNumberDistribution',char(Mutants{j}),'.eps'));

    Var2{j} = strcat(num2str(X(1),2)); 
    Var3{j} = strcat(num2str(X(2),2)); 
    Var4{j} = strcat(num2str(X(3),2)); 
    Var5{j} = strcat(num2str(X(4),2)); 
    Var6{j} = strcat(num2str(X(5),2)); 
    Var7{j} = strcat(num2str(X(6),2)); 
    Var8{j} = strcat(num2str(X(7),2)); 
    
    
    %Plot number histogram
    edges = (0:10)-0.5;
    figure(400+j)
    hold on
    title(strcat('Nr. short-lived bulbous  tips (',char(Mutants{j}),')'),'FontSize',14)
    [counts,edges] = histcounts(data_s,edges);
    X = counts./sum(counts);
    histogram('BinEdges',edges,'BinCounts',X)
    line([m_s  m_s],[0 0.2],'Color','r','LineStyle','--','LineWidth',3)
    line([m_s+s_s m_s+s_s],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    line([m_s-s_s  m_s-s_s],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    ylim([0 0.9])
    ylabel('Frequency/Probability','FontWeight','bold','FontSize',14)
    xlabel('Number bulbous tips')
    %print(400+j,'-depsc2',strcat('./Figures/ShortBulbousNumberDistribution',char(Mutants{j}),'.eps'));

    Var2_2{j} = strcat(num2str(X(1),2)); 
    Var2_3{j} = strcat(num2str(X(2),2)); 
    Var2_4{j} = strcat(num2str(X(3),2)); 
    Var2_5{j} = strcat(num2str(X(4),2)); 
    Var2_6{j} = strcat(num2str(X(5),2)); 
    Var2_7{j} = strcat(num2str(X(6),2)); 
    Var2_8{j} = strcat(num2str(X(7),2)); 
    
end
        
Names = {'Mutant';'Nr_syn_Bulbous'};
%print data to table
T = table(Mutants,Var1);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Number (standard deviation) of synaptogenic bulbous tips per growth cone (P60) and time instance';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);    

Names = {'Mutant';'Nr_s_Bulbous'};
%print data to table
T = table(Mutants,Var2);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Number (standard deviation) of short-lived bulbous tips per growth cone (P60) and time instance';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);    


%%

Names = {'Mutant';'n0';'n1';'n2';'n3';'n4';'n5';'n6'};
%print data to table
T = table(Mutants,Var2,Var3,Var4,Var5,Var6,Var7,Var8);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Probability the the number of synaptogenic bulbous tips per growth cone (P60) and time instance is n';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);  

Names = {'Mutant';'n0';'n1';'n2';'n3';'n4';'n5';'n6'};
%print data to table
T = table(Mutants,Var2_2,Var2_3,Var2_4,Var2_5,Var2_6,Var2_7,Var2_8);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Probability the the number of transient bulbous tips per growth cone (P60) and time instance is n';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);  


