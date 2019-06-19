function [] = AnalyzeBulbousLifeTimes()

% This routine depicts the lifetime of bulbous tips and computes parameter
% c4 'death rate parameter' of bulbous tips

Mutants = {'Mutant';'WT';'DLar';'LiprinA';'Syd1';'Trio'};
LifeTimeThreshold = 40;

%Bulbous tip life times (min)
Data.BulbTip.WT = [60 28 53 4 58 8 60 1 60 60 40 60 56];
Data.BulbTip.DLar = [3 1 9 4 2 60 1];
Data.BulbTip.LiprinA = [6 2 8 9 1 6 43 32 1 2 6 8 7 29 9];
Data.BulbTip.Syd1 = [2 7 3 4 41 1 1 9 10 3 12 9 1 1 1 2 12 23 3 1 6 17 2 1 3 1 1 2 59 10 27 2 5 20 15 2 1 1 10 1 1 10 60 13 4 8 2 1 2];
Data.BulbTip.Trio = [4 59 10 4 5 15 2 8 4 1 3 60 12 18 6 2 6 53 1 2 11 3 15 3 2 25 43 21];

%corresponding growth cones
Data.GC.WT =       [1   1   2   3  3   3   4   4   5   5   6   6   6];%growth cone
Data.GC.DLar =         [1 1   2   3   3   3   3];
Data.GC.LiprinA =       [1  1  1  1  1  1  2  2  2  2  3  3  3  3   3];
Data.GC.Syd1 =       [1 1  1  1  1  1 1 1  1  1 1  1  1   2  2  2  3  3  3  3  3   3 3  3 3  4 4 4  4 4  4  4 4 4  4  4 4 4 4  4 4 5  5   5 5  5  5 5 5];
Data.GC.Trio =       [1 1   1 1     1   1 1 1   1  2  2  2  2  2  2 3 3  3  3 4 4   4 4  4  4   4  4 4];

edges = (0:10:60);
Labels = cellstr(num2str(edges'));
Labels{end} = strcat('>',Labels{end});
Colors = ['g';'r';'k';'m';'y';'c'];

for i = 2:length(Mutants)
    mutant = char(Mutants(i))
    BulbTips = Data.BulbTip.(mutant);
    GrowthCones = Data.GC.(mutant);
    
    %make a histogram of the life times
    figure(i-1)
    [N,edges] = histcounts(BulbTips,edges);%
    heights = N./sum(N);
    xplot = edges(1:end-1)+0.5*(edges(2)-edges(1));
    bar(xplot,heights,1)
    title(mutant)
    xlabel('Lifetime (min)')
    ylabel('Frequency')
    fs = 18;
    set(gca,'FontSize',fs);
    set(get(gca,'title'),'Fontsize',fs);  
    set(get(gca,'xlabel'),'Fontsize',fs); 
    set(get(gca,'ylabel'),'Fontsize',fs);   
    set(gca,'XTick',edges,'XTickLabels',Labels)
    xlim([0 60])

    NrGC = length(unique(GrowthCones));
    NrInCategory = nan(length(N),NrGC);
    
    %plot the individual numbers per life time category and growth cone
    for j = 1:NrGC
        idx = find(GrowthCones==j);
        NrBulbsPerGC = length(idx);
        
        %measurement of short lived bulbs, assuming all disappear
        %eventually
        BulbsInGC = BulbTips(idx);

        for z = 1:length(N)-1
            NrInCategory(z,j) = sum(BulbsInGC>=edges(z) & BulbsInGC<edges(z+1));
        end
        z = length(N);
        NrInCategory(z,j) = sum(BulbsInGC>=edges(z) & BulbsInGC<=edges(z+1));
        txt = strcat(num2str(NrInCategory(:,j)),'/',num2str(NrBulbsPerGC));
        for z = 1:length(N)
           text(xplot(z),0.5-0.1*(j-1) ,txt(z,:),'Color',Colors(j),'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',14)
        end
    end
  
    %print(i-1,'-dpdf',strcat(Path,Name,'.pdf'))
    %print(i-1,'-dtiff',strcat(Path,Name,'.tiff'))
    %print(i-1,'-depsc',strcat(Path,Name,'.eps'))
    
    lgx = BulbTips>=LifeTimeThreshold;
    if mean(BulbTips(lgx)) > mean(BulbTips(~lgx))
        Long = lgx;
        Short = ~lgx;
    else
        Long = ~lgx;
        Short = lgx;
    end
    
    ShortLived = BulbTips(Short);
    
    if i > 2
        c4 = 1./mean(ShortLived)
    else
        c4 = 1/120
    end
    MeanLifetime = mean(ShortLived)
    StdLiveTime = std(ShortLived)

end
