function [Data]= DataImport()
% Imports the tabularized data regarding the live-tracking of each
% filopodium at P40 and P60 respectively for the different genotypes and
% saves it to ./FiloData/AllData.mat. It also generates a number of plots
% to visualize the live-tracking data.
% The data is imported into a Data structure of the form Data.<mutant>.<time>.<type>
% mutant: 'DLar' or 'LiprinA' or 'Syd1' or 'Trio' or 'WT'
% time: 'P40' or 'P60'
% type: 'sF' or 'ellF'  (short- and long-lived filopodia)
% Example: Data.DLar.P40.sF contains all data of short-lived filopodia at
% P40 for the DLar knockout mutant.
%The following data was stored 'LTimes' -Lifetimes; 'StartTimes' -when the 
% filopodium emerged; 'EndTimes' -when it disappeared; 'GC' - which growth
% cone it belongs to 

% 1: Import data from the mutants
path = './FiloData/';
Mutants = {'DLar'; 'LiprinA'; 'Syd1';'Trio';'WT'};

LClassificationThreshold = 8% in min shortlived have lifetime smaller than this threshold  change in lines 113 and 464

for j = 1:length(Mutants)
    mutant = char(Mutants(j));
    idx2 = [0 0];
    %Path name and entering the path
    FullPath = strcat(path,mutant);
    cd(FullPath)
    %check content and find files
    k = ls;
    if j ~= 5 % Mutants
       FileTypePos = strfind(k,'.xlsx');
       incr1 = 4;
       incr2 = 6;
       gcID = 'GC';
       timeID = 'fast3';
       fieldnameStart = 'StartTime';
       fieldnameEnd = 'EndTime';
       fieldnameLT = 'LifeTime';
    else % WT
       FileTypePos = strfind(k,'.csv'); 
       incr1 = 3;
       incr2 = 5;
       gcID = 'gc';
       timeID = 'P60';
       fieldnameStart = 'StartTimeStep';
       fieldnameEnd = 'EndTimeStep';
       fieldnameLT = 'LifeTime_min_';
    end
 
    NrFiles = length(FileTypePos);
    Filenames = cell(NrFiles,1); % initialize empty cell string 
    %Initialize growth cone assignment
    GCs = nan(NrFiles,1);
    
    pos = 1;
    % store filenames in repository
    for i = 1:NrFiles
        filename = k(pos:FileTypePos(i)+incr1);
        Filenames{i} = filename;
        pos = FileTypePos(i)+incr2;
        
        %assign growth cone
        GC_pos = strfind(char(filename),gcID)+2;
        GCs(i) = str2num(filename(GC_pos));
    end
        %Asssign files to P40 or P60
        Assignment = zeros(1,NrFiles);
        for z = 1:length(Filenames)
            filename = char(Filenames{z});
            if contains(filename,timeID)
                %P60
                Assignment(z) = 2;
            else
                %P40
                Assignment(z) = 1;
            end
        end

        % save statistics for each mutant
        NrP40 = sum(Assignment == 1); % Nr of P40 files
        NrP60 = sum(Assignment == 2); % Nr of P60 files
        FiloInfoP40 = nan(300,4,NrP40);%For each Filo (dim1): Start, End, Lifetime (dim2), dim3 = growth cone
        FiloInfoP40_sF = nan(300,4,NrP40);
        FiloInfoP40_ellF = nan(300,4,NrP40);
        FiloInfoP60 = nan(300,4,NrP60);%For each Filo (dim1): Start, End, Lifetime (dim2), dim3 = growth cone
        FiloInfoP60_sF = nan(300,4,NrP60);
        FiloInfoP60_ellF = nan(300,4,NrP60);

        LifeTimeValues_sF = nan(500,max([NrP40 NrP60]),2);
        LifeTimeValues_ellF = nan(500,max([NrP40 NrP60]),2);
        
        NbrsFilosAllExperimentsP40_ellF = zeros(NrP40,60);
        NbrsFilosAllExperimentsP40_sF = zeros(NrP40,60);
        NbrsFilosAllExperimentsP60_ellF = zeros(NrP60,60);
        NbrsFilosAllExperimentsP60_sF = zeros(NrP60,60);
        
        counterP40 = 0;
        counterP60 = 0;
        for z = 1:length(Filenames)
            filename = char(Filenames{z});
            %read table
            T = readtable(filename);
            NrEntries = length(T.(fieldnameStart));

            StartTimes = T.(fieldnameStart);
            EndTimes = T.(fieldnameEnd);
            LifeTimes = T.(fieldnameLT);

            %save the data
            idx1 = Assignment(z); % p40 or P60
            idx2(idx1) = idx2(idx1)+1;
            idx3 = idx2(idx1);%indicates the xth file of either P40 or P60

            TF_fast = (LifeTimes < LClassificationThreshold);
            TF_stable = (LifeTimes >= LClassificationThreshold);
            
            if idx1 == 1
                FiloInfoP40(1:NrEntries,1,idx3) = StartTimes;
                FiloInfoP40(1:NrEntries,2,idx3) = EndTimes;
                FiloInfoP40(1:NrEntries,3,idx3) = LifeTimes;
                FiloInfoP40(1:NrEntries,4,idx3) = ones(NrEntries,1).*GCs(z);
                
                FiloInfoTmp1 = FiloInfoP40(TF_fast,:,idx3);
                FiloInfoTmp2 = FiloInfoP40(TF_stable,:,idx3);
            
                FiloInfoP40_sF(1:sum(TF_fast),:,idx3) = FiloInfoTmp1;
                FiloInfoP40_ellF(1:sum(TF_stable),:,idx3) = FiloInfoTmp2;
                
                
            else
                FiloInfoP60(1:NrEntries,1,idx3) = StartTimes;
                FiloInfoP60(1:NrEntries,2,idx3) = EndTimes;
                FiloInfoP60(1:NrEntries,3,idx3) = LifeTimes;
                FiloInfoP60(1:NrEntries,4,idx3) = ones(NrEntries,1).*GCs(z);
                
                FiloInfoTmp1 = FiloInfoP60(TF_fast,:,idx3);
                FiloInfoTmp2 = FiloInfoP60(TF_stable,:,idx3);
            
                FiloInfoP60_sF(1:sum(TF_fast),:,idx3) = FiloInfoTmp1;
                FiloInfoP60_ellF(1:sum(TF_stable),:,idx3) = FiloInfoTmp2;
                
            end
            

            LifeTimesVec = LifeTimes';

            % Classify filopodia into different types based on their
            % lifetime
            TF = LifeTimesVec >= LClassificationThreshold;
            LifeTimes_sF= LifeTimesVec(~TF);
            LifeTimes_ellF = LifeTimesVec(TF);
            LifeTimeValues_sF(1:length(LifeTimes_sF),idx3,idx1) = LifeTimes_sF';
            LifeTimeValues_ellF(1:length(LifeTimes_ellF),idx3,idx1) = LifeTimes_ellF';

             %% Compute a number of statistics
            %Plot Filo & Bulbous Numbers

            if idx1 == 1
                TmpLGX = ~isnan(FiloInfoP40(:,1,idx3));%trim vector
                FiloInfo = FiloInfoP40(TmpLGX,:,idx3); 
            else
                TmpLGX = ~isnan(FiloInfoP60(:,1,idx3));%trim vector
                FiloInfo = FiloInfoP60(TmpLGX,:,idx3);
            end
            %Only real measurements (arose and disappeared)
            TmpLGX = ~isinf(FiloInfo(:,3));
            FiloInfo = FiloInfo(TmpLGX,:); 
            %Assign what is filopodium and what is stab. Filo. 
            FiloLGX_sF = FiloInfo(:,3) < LClassificationThreshold;
            FiloLGX_ellF = FiloInfo(:,3) >= LClassificationThreshold;
            FiloInfoTmp_sF = FiloInfo(FiloLGX_sF,:);
            FiloInfoTmp_ellF = FiloInfo(FiloLGX_ellF,:);

            NrFilos_sF = sum(FiloLGX_sF);
            NrFilos_ellF = sum(FiloLGX_ellF);
            NrFilosPerMinute_sF = zeros(1,60);
            for counter = 1:NrFilos_sF
                Start = FiloInfoTmp_sF(counter,1)+1;
                End = FiloInfoTmp_sF(counter,2)+1;
                NrFilosPerMinute_sF(Start:End) = NrFilosPerMinute_sF(Start:End)+1;
            end
            NrFilosPerMinute_ellF = zeros(1,60);   
            for counter = 1:NrFilos_ellF
                Start = FiloInfoTmp_ellF(counter,1)+1;
                End = FiloInfoTmp_ellF(counter,2)+1;
                NrFilosPerMinute_ellF(Start:End)= NrFilosPerMinute_ellF(Start:End)+ 1;
            end
            figure(101+j)
            subplot(2,3,(idx1-1)*3+1)
            hold on
            plot(NrFilosPerMinute_sF,':','LineWidth',2)
            subplot(2,3,(idx1-1)*3+2)
            hold on
            plot(NrFilosPerMinute_ellF,':','LineWidth',2)
            %------------------
            
            if idx1 == 1
                counterP40 = counterP40 +1;
                NbrsFilosAllExperimentsP40_ellF(counterP40,:) = NrFilosPerMinute_ellF;
                NbrsFilosAllExperimentsP40_sF(counterP40,:) = NrFilosPerMinute_sF;
            else
                counterP60 = counterP60 +1;
                NbrsFilosAllExperimentsP60_ellF(counterP60,:) = NrFilosPerMinute_ellF;
                NbrsFilosAllExperimentsP60_sF(counterP60,:) = NrFilosPerMinute_sF;
            end

            
        end %end over files
        
        %save the data to a structure 
        % P40, sF
        FiloInfo_tmp = FiloInfoP40_sF(:,:,1);
        if NrP40 > 1 
            for i = 2:NrP40
                FiloInfo_tmp = cat(1,FiloInfo_tmp,FiloInfoP40_sF(:,:,i));
            end
        end
        LGX = ~isnan(FiloInfo_tmp(:,1));
        FiloInfo_tmp2 = FiloInfo_tmp(LGX,:);
        Data.(mutant).P40.sF.LTimes = FiloInfo_tmp2(:,3);    
        Data.(mutant).P40.sF.StartTimes = FiloInfo_tmp2(:,1); 
        Data.(mutant).P40.sF.EndTimes = FiloInfo_tmp2(:,2); 
        Data.(mutant).P40.sF.GC = FiloInfo_tmp2(:,4); 
        % P60, sF
        FiloInfo_tmp = FiloInfoP60_sF(:,:,1);
        if NrP40 > 1 
            for i = 2:NrP40
                FiloInfo_tmp = cat(1,FiloInfo_tmp,FiloInfoP60_sF(:,:,i));
            end
        end
        LGX = ~isnan(FiloInfo_tmp(:,1));
        FiloInfo_tmp2 = FiloInfo_tmp(LGX,:);
        Data.(mutant).P60.sF.LTimes = FiloInfo_tmp2(:,3);    
        Data.(mutant).P60.sF.StartTimes = FiloInfo_tmp2(:,1); 
        Data.(mutant).P60.sF.EndTimes = FiloInfo_tmp2(:,2); 
        Data.(mutant).P60.sF.GC = FiloInfo_tmp2(:,4); 
        % P40, ellF
        FiloInfo_tmp = FiloInfoP40_ellF(:,:,1);
        if NrP40 > 1 
            for i = 2:NrP40
                FiloInfo_tmp = cat(1,FiloInfo_tmp,FiloInfoP40_ellF(:,:,i));
            end
        end
        LGX = ~isnan(FiloInfo_tmp(:,1));
        FiloInfo_tmp2 = FiloInfo_tmp(LGX,:);
        Data.(mutant).P40.ellF.LTimes = FiloInfo_tmp2(:,3);    
        Data.(mutant).P40.ellF.StartTimes = FiloInfo_tmp2(:,1); 
        Data.(mutant).P40.ellF.EndTimes = FiloInfo_tmp2(:,2); 
        Data.(mutant).P40.ellF.GC = FiloInfo_tmp2(:,4); 
        % P60, sF
        FiloInfo_tmp = FiloInfoP60_ellF(:,:,1);
        if NrP40 > 1 
            for i = 2:NrP40
                FiloInfo_tmp = cat(1,FiloInfo_tmp,FiloInfoP60_ellF(:,:,i));
            end
        end
        LGX = ~isnan(FiloInfo_tmp(:,1));
        FiloInfo_tmp2 = FiloInfo_tmp(LGX,:);
        Data.(mutant).P60.ellF.LTimes = FiloInfo_tmp2(:,3);    
        Data.(mutant).P60.ellF.StartTimes = FiloInfo_tmp2(:,1); 
        Data.(mutant).P60.ellF.EndTimes = FiloInfo_tmp2(:,2); 
        Data.(mutant).P60.ellF.GC = FiloInfo_tmp2(:,4); 
        
        
    cd ..
    cd ..

    %histogram of the number distribution over the 60 min live-imaging interval
    figure(201+j)
    hold on
    edges = (0:20)-0.5;
    for i = 1:4
        subplot(2,2,i)
         hold on
        switch i
            case 1
                Dta = NbrsFilosAllExperimentsP40_sF(:);
               title(strcat('Nr. instab filo.(',mutant,')'),'FontSize',16)
               ylabel('P40','FontSize',16)
            case 2
                 title(strcat('Nr. stab. filo.(',mutant,')'),'FontSize',16)
                 Dta = NbrsFilosAllExperimentsP40_ellF(:);
            case 3
               ylabel('P60','FontSize',16)
               Dta = NbrsFilosAllExperimentsP60_sF(:);
               xlabel('number','FontSize',14)
            case 4
                Dta = NbrsFilosAllExperimentsP60_ellF(:);
               xlabel('number','FontSize',14)
        end
         [counts,edges] = histcounts(Dta,edges);
        histogram('BinEdges',edges,'BinCounts',counts./sum(counts))
        MeanL = mean(Dta);
        StdL = std(Dta); 
        line([MeanL  MeanL],[0 0.2],'Color','r','LineStyle','--','LineWidth',3)
        line([MeanL+StdL MeanL+StdL],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
        line([MeanL-StdL  MeanL-StdL],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
        ylim([0 0.3])
        
    end
   
    %print(102+j,'-dpng',strcat('FiloNumberDistributions',mutant,'.png'))
    
    figure(101+j)
    subplot(2,3,1)
    times = 1:60;
    %Nr Filo P40
    Means1 = mean(NbrsFilosAllExperimentsP40_sF);
    TF = isoutlier(Means1);
    Means1(TF) = nan;
    plot(times,Means1,'k-','LineWidth',4)
    ylim([0 15])
    %Nr Filo P60
    Means2 = mean(NbrsFilosAllExperimentsP60_sF);
    TF = isoutlier(Means2);
    Means2(TF) = nan;
    subplot(2,3,4)
    plot(times,Means2,'k-','LineWidth',4)
    ylim([0 15])
    %Nr Stab. Filo P40
    Means3 = mean(NbrsFilosAllExperimentsP40_ellF);
    TF = isoutlier(Means3);
    Means3(TF) = nan;
    subplot(2,3,2)
    plot(times,Means3,'k-','LineWidth',4)
    ylim([0 20])
    %Nr Stab. Filo P60
    Means4 = mean(NbrsFilosAllExperimentsP60_ellF);
    TF = isoutlier(Means4);
    Means4(TF) = nan;
    subplot(2,3,5)
    plot(times,Means4,'k-','LineWidth',4)
    ylim([0 20])
    %Ratio in P40
    Ratio1 = Means3./Means1;% correction for non-detection
    TF = isoutlier(Ratio1,'movmedian',5);
    Ratio1(TF) = nan;
    Ratio1(isnan(Ratio1)) = inf;
    TF2 = isinf(Ratio1);
    Ratio1 = Ratio1(~TF2);
    subplot(2,3,3)
    plot(times(~TF2),Ratio1,'k-','LineWidth',4)
    hold on
    subplot(2,3,3)
    Ratio1(isnan(Ratio1)) = inf;
    Ratio1 = Ratio1(~isinf(Ratio1));
    line([1 60],[mean(Ratio1) mean(Ratio1)],'LineStyle',':','Color','r','LineWidth',3)
    ylim([0 10])
    title('Ratio: stab.Filo-to-instab.Filo')

    %Ratio in P60
    Ratio2 = Means4./Means2;% correction for non-detection
    TF1 = isoutlier(Ratio2,'movmedian',5);
    Ratio2(TF1) = nan;
    Ratio2(isnan(Ratio2)) = inf;
    TF2 = isinf(Ratio2);
    Ratio2 = Ratio2(~TF2);
    subplot(2,3,6)
    plot(times(~TF2),Ratio2,'k-','LineWidth',4)
    line([1 60],[mean(Ratio2) mean(Ratio2)],'LineStyle',':','Color','r','LineWidth',3)
    ylim([0 10])
    %title('Ratio: stab.Filo-to-Filo')

    figure(101+j)
    subplot(2,3,1)
    ylabel('P40','FontWeight','bold','FontSize',14)
    title(strcat('Nr. instab. filo. (',mutant,')'))
    subplot(2,3,4)
    ylabel('P60','FontWeight','bold','FontSize',14)
    subplot(2,3,2)
    title(strcat('Nr. stab. filo.(',mutant,')'))
    subplot(2,3,4)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    subplot(2,3,5)
    xlabel('time  (min)','FontWeight','bold','FontSize',12)
    subplot(2,3,6)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    %print(101+j,'-dpng',strcat('FiloNumbers',mutant,'.png'))

end % end over Mutants
try
    delete(strcat(path,'AllData.mat'));
catch
end
save(strcat(path,'AllData.mat'),'Data');