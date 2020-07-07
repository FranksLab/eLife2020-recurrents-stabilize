clearvars
close all
clc

%% Parameters
%% PCX original data
Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_bulb_awk_kx_F.txt';
moff = 1;
normtodiff = 0;

clear ndcca pcca

for rewi = 1:2
    
    BinSizes = [0.5];
    if rewi == 1
        PST = [0 0.5];
    else
        PST = [-2 -1.5];
    end
    statelist = ['a'];
    
    % Do celltyping.
    specialparams.FRlim = 1/100;
    specialparams.UVlim = 50;
    specialparams.DFRlim = 100;
    
    [TypeIdxa, TypeStack] = CellTyper (Catalog, 'Stable', specialparams);
    % Load in all EFD data to save on processing time.
    efds = EFDloader(Catalog);
    % Concatenate the catalogued data [X,Y,n_bins] = PseudoPopulator(Catalog, BinSizes, PST)
    [Xa,Ya,n_binsa,trda] = CrossPopulatorFifteenTrials(Catalog, efds, BinSizes, PST, statelist); % if abs(BinSize-diff(PST))<BinSize there will be 1 bin.
    
    
    rset{1} = [4,6,7,9];
    rset{2} = [5,8,10,12];
    rset{3} = [11];
    
    %%
    VOI = 2:7;
    
    for r = 1:length(rset)
        Xan = trda(rset{r},:);
        
        for k = 1:size(Xan,1)
            Xan{k,1} = Xan{k,1}(:,TypeIdxa{rset{r}(k),1});
            
            YnVOI = Ya{1}(ismember(Ya{1},VOI));
            
            nd = 3;
            Xkz = Xan{k,1}(ismember(Ya{1},VOI),:);
            Xkz = zscore(Xkz);
            [mappedA,mapping] = compute_mapping(Xkz,'PCA',nd);
            
            
            n = squareform(pdist(Xan{k,1}(ismember(Ya{1},VOI),:)));
            for tr = 1:15
                interoa{r,rewi}(tr,k) = nanmean(reshape(n(tr:15:end,tr:15:end),[],1));
            end
            
            for v = 1:length(VOI)
                
                %% PCA set for illustration
                % PC distance
                currset = mappedA(ismember(YnVOI,VOI(v)),:);
                
                if rset{r}(k) == 10 && rewi ==1
                    size(Xkz)
                    %%
                    
                    figure(9)
                    clf
                    printpos([100 00 800 550])
                    
                    figure(10)
                    clf
                    printpos([100 00 900 550])
                    
                    littlevlist = [2,3,6];
                    for vl = 1:length(littlevlist)
                        currset = mappedA(ismember(YnVOI,VOI(littlevlist(vl))),:);
                        %
                        CM = circshift(copper(32),vl,2);
                        CM = CM(10:end,:);
                        %
                        
                        figure(10)
                        ax = subplot(1,3,0+vl);
                        colormap(ax,CM)
                        covarianceEllipse3D(mean(currset(11:15,[1,2,3])),cov(currset(11:15,[1,2,3])),15)
                        hold on
                        plot3(currset(:,1),currset(:,2),currset(:,3),'--','color',CM(10,:))
                        scatter3(currset(1:15,1),currset(1:15,2),currset(1:15,3),250,1:15,'.')
                        
                        xlabel('PC1')
                        ylabel('PC2')
                        zlabel('PC3')
                        xlim([-4 8]); ylim([-6 6]); zlim([-4 4]);
                        set(gca,'XDir','reverse','clipping','off')
                        axis square
                        view([-15 25])
                    end
                    %%
                end
                
                % neural distance
                currset = Xan{k,1}(ismember(Ya{1},VOI(v)),:);
                ratesa{r,rewi}(:,k,v) = nanmean(currset,2);
                
                currset_trunc = currset(1:15,:);
                currset_std = (currset(11:15,:));
                
                te = pdist2(currset_trunc, currset_std, 'euclidean'); % distances from standard
                te(te ==0) = nan;
                
                ndccal{r,rewi}(:,k,v) = nanmean(te,2); % distances from standard
                
                temp = squareform(pdist(currset,'euclidean'));
                llt = 11:15;
                tlt = tril(temp(llt,llt));
                tlt(tlt == 0) = nan;
                ltvara{r,rewi}(:,k,v) = nanmean(tlt(:));
                
                
                % sniffcount
                T = readtable(Catalog, 'Delimiter', ' ');
                
                if strcmp(T.VOI(rset{r}(k)),'A')
                    VOIe = [4,7,8,12,15,16];
                elseif strcmp(T.VOI(rset{r}(k)),'C')
                    VOIe = [11,7,8,6,12,10];
                end
                
                snifftimes = efds(rset{r}(k)).ValveTimes.PREXTimes{VOIe(v)}(1:15);
                for tr  = 1:15
                    sniffidxs = find(efds(rset{r}(k)).PREX >= (snifftimes(tr) + PST(1)) & efds(rset{r}(k)).PREX <= snifftimes(tr) + PST(2),1);
                    if isempty(sniffidxs)
                        sniffidxs = find(efds(rset{r}(k)).PREX >= (snifftimes(tr) + PST(1)-1) & efds(rset{r}(k)).PREX <= snifftimes(tr) + PST(2),1);
                    end
                    sniffsa{r,rewi}(tr,k,v) = nanmedian(1./efds(rset{r}(k)).BreathStats.BbyB.Width(sniffidxs));
                end
                
            end
        end
    end
end

%% plot curves
figure(68)
clf
printpos([300 200 700 500])
colores = [.5 .4 1; 1 .4 .4];
r = 1;
clear combinedneuraldistance combinedsniffs
for rewi = 2:-1:1
    subplot(1,2,1)
    combinedneuraldistance{r,rewi} = bsxfun(@minus,cat(2,[ndccal{:,rewi}]), cat(2,[ltvara{:,rewi}])) ;% ./ cat(2,[ndccal{r,2},ndccnl{r,2}]);
    
    if rewi == 1
        hold on
        plot(mean(combinedneuraldistance{r,rewi}(:,:),2),'o','markersize',4,'markerfacecolor','k','markeredgecolor','none');
        errorbar(1:15,mean(combinedneuraldistance{r,rewi}(:,:),2),sem(combinedneuraldistance{r,rewi}(:,:)'),'.','color','k')
    else
        boundedline(1:15,mean(combinedneuraldistance{r,rewi}(:,:),2),sem(combinedneuraldistance{r,rewi}(:,:)'),'k','cmap',[0 0 0])
    end
    hold on
    axis square
    box off
    ylim([-1 6])
    ylabel('norm. neural distance from standard')
    xlim([0 16])
    set(gca,'XTick',1:2:15)
    xlabel('trials')
    %
    subplot(1,2,2)
    yyaxis left
    combinedsniffs{r,rewi} = cat(2,[sniffsa{:,rewi}]);
    
    
    if rewi == 1
        hold on
        plot(mean(combinedsniffs{r,rewi}(:,:),2),'o','markersize',4,'markerfacecolor',colores(1,:),'markeredgecolor','none');
        errorbar(1:15,mean(combinedsniffs{r,rewi}(:,:),2),sem(combinedsniffs{r,rewi}(:,:)'),'.','color',colores(1,:))
    else
        boundedline(1:15,mean(combinedsniffs{r,rewi}(:,:),2),sem(combinedsniffs{r,rewi}(:,:)'),'k','cmap',colores(1,:))
    end
    hold on
    axis square
    box off
    xlim([0 16])
    ylim([1.5 5])
    ylabel('sniff rate during odor')
    set(gca,'XTick',1:2:15,'YTick',2:5)
    set(gca,'clipping','off')
    
    yyaxis right
    combinedneuraldistance{r,rewi} = cat(2,[ratesa{:,rewi}])/.5;
    if rewi == 1
        hold on
        plot(mean(combinedneuraldistance{r,rewi}(:,:),2),'o','markersize',4,'markerfacecolor',colores(2,:),'markeredgecolor','none');
        errorbar(1:15,mean(combinedneuraldistance{r,rewi}(:,:),2),sem(combinedneuraldistance{r,rewi}(:,:)'),'.','color',colores(2,:))
    else
        boundedline(1:15,mean(combinedneuraldistance{r,rewi}(:,:),2),sem(combinedneuraldistance{r,rewi}(:,:)'),'k','cmap',colores(2,:))
    end
    hold on
    axis square
    box off
    xlim([0 16])
    ylim([2.5 10])
    ylabel('mean pop. firing rate')
    set(gca,'clipping','off')
    set(gca,'XTick',1:2:15)
    xlabel('trials')
    
end

%% regression
figure(3423)
printpos([300 300 600 300])
colores = parula(4);
clf

for rewi = 1%:2
    
    
    % for distance to last trial
    combinedneuraldistance = ((cat(2,[ndccal{:,rewi}])));
    
    trialnum = repmat([1:15]',1,size(combinedneuraldistance,2),6);
    
    
    combinedneuralrates = cat(2,[ratesa{:,rewi}]);
    combinedneuralrates = combinedneuralrates(1:15,:,:);
    %     combinedneuralrates = (bsxfun(@minus, combinedneuralrates(1:end-1,:,:), combinedneuralrates(end,:,:)));
    
    combinedsniffs = cat(2,[sniffsa{:,rewi}]);
    combinedsniffs = combinedsniffs(1:15,:,:);
    
    
    mdl = fitlm([zscore(trialnum(:)), zscore(combinedneuralrates(:)), zscore(combinedsniffs(:))],(combinedneuraldistance(:)), ...
        'VarNames',{'trials','rates','sniffs','distance'});
    
    subplot(1,2,1)
    if rewi == 1
        bar(1:length(mdl.Coefficients.Estimate(2:end)),mdl.Coefficients.Estimate(2:end),.8,'b')
        
    end
    hold on
    errorbar(1:length(mdl.Coefficients.Estimate(2:end)),mdl.Coefficients.Estimate(2:end),mdl.Coefficients.SE(2:end),'k.')
    set(gca,'XTick',1:length(mdl.Coefficients.Estimate(2:end)),'XTickLabel',mdl.CoefficientNames(2:end))
    box off
    xlim([0 length(mdl.Coefficients.Estimate(2:end))+1])
    ax = gca;
    ax.XTickLabelRotation = 60;
    hold on
    ylabel('regression coefficient')
    ylim([-1 3])
    
    %     Residuals plot (control for sniff and rate)
    mdl = fitlm([zscore(combinedneuralrates(:)), zscore(combinedsniffs(:))],combinedneuraldistance(:), ...
        'VarNames',{'rates','sniffs','distance'});
    distanceresiduals = mdl.Residuals.Raw;
    distanceresiduals = reshape(distanceresiduals, size(combinedneuraldistance));
    %
    for r = 1
        subplot(1,2,2)
        dr = distanceresiduals(:,:,:);
        %
        if rewi == 1
            hold on
            errorbar(1:15,mean(dr(:,:),2),sem(dr(:,:)'),'.','color',colores(r,:))
            plot(mean(dr(:,:),2),'o','markersize',4,'markerfacecolor',colores(r,:),'markeredgecolor','none');
        else
            errorbar(1:15,mean(dr(:,:),2),sem(dr(:,:)'),'.','color',colores(r,:))
            plot(mean(dr(:,:),2),'o','markersize',4,'markeredgecolor',colores(r,:),'markerfacecolor','w');
        end
        box off
        axis square
        xlim([0 16])
    end
    
    
end



