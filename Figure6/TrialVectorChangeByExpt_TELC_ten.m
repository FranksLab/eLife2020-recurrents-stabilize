clearvars
close all
clc

%% Parameters
%% PCX original data
Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_pcx_awk_kx_telc_stabletrials.txt';
rset{1} = [2,4,5,6];

% Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_pcx_awk_kx_telcontrol_stabletrials.txt';
% rset{1} = [1:4];

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
    [Xa,Ya,n_binsa,trda] = CrossPopulatorTenTrials(Catalog, efds, BinSizes, PST, statelist); % if abs(BinSize-diff(PST))<BinSize there will be 1 bin.
    
    %%
    VOI = 2:7;
    
    for r = 1:length(rset)
        Xan = trda(rset{r},:);
        
        for k = 1:size(Xan,1)
            Xan{k,1} = Xan{k,1}(:,TypeIdxa{rset{r}(k),1});
            YnVOI = Ya{1}(ismember(Ya{1},VOI));
            n = squareform(pdist(Xan{k,1}(ismember(Ya{1},VOI),:)));
            
            for v = 1:length(VOI)
                
                % neural distance
                currset = Xan{k,1}(ismember(Ya{1},VOI(v)),:);
                ratesa{r,rewi}(:,k,v) = nanmean(currset,2);
                
                currset_trunc = currset(1:10,:);
                currset_std = (currset(7:10,:));
                
                te = pdist2(currset_trunc, currset_std, 'euclidean'); % distances from standard
                te(te ==0) = nan;
                
                ndccal{r,rewi}(:,k,v) = nanmean(te,2); % distances from standard
                
                temp = squareform(pdist(currset,'euclidean'));
                llt = 7:10;
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
                
                snifftimes = efds(rset{r}(k)).ValveTimes.PREXTimes{VOIe(v)}(1:10);
                for tr  = 1:10
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

%% regression
figure(3423)
printpos([300 300 600 300])
colores = parula(4);
clf

rewi = 1;
trialnum = repmat([1:10]',1,length(rset{1}),6);

combinedneuraldistance = ((cat(2,[ndccal{:,rewi}])));

combinedneuralrates = cat(2,[ratesa{:,rewi}]);
combinedneuralrates = combinedneuralrates(1:10,:,:);

combinedsniffs = cat(2,[sniffsa{:,rewi}]);
combinedsniffs = combinedsniffs(1:10,:,:);

mdl = fitlm([zscore(trialnum(:)), zscore(combinedneuralrates(:)), zscore(combinedsniffs(:))],combinedneuraldistance(:), ...
    'VarNames',{'trials','rates','sniffs','distance'});

subplot(1,2,rewi)
bar(1:length(mdl.Coefficients.Estimate(2:end)),mdl.Coefficients.Estimate(2:end),.8,'b')
hold on
plot(1:length(mdl.Coefficients.pValue(2:end)),(mdl.Coefficients.pValue(2:end)<.05)*8,'*')
ylim([-4 8])
hold on
errorbar(1:length(mdl.Coefficients.Estimate(2:end)),mdl.Coefficients.Estimate(2:end),mdl.Coefficients.SE(2:end),'k.')
set(gca,'XTick',1:length(mdl.Coefficients.Estimate(2:end)),'XTickLabel',mdl.CoefficientNames(2:end))
box off
xlim([0 length(mdl.Coefficients.Estimate(2:end))+1])
ax = gca;
ax.XTickLabelRotation = 60;
hold on
ylabel('regression coefficient')

mdl = fitlm([zscore(combinedneuralrates(:)), zscore(combinedsniffs(:))],combinedneuraldistance(:), ...
    'VarNames',{'rates','sniffs','distance'});
distanceresiduals = mdl.Residuals.Raw;
distanceresiduals = reshape(distanceresiduals, size(combinedneuraldistance));

r = 1;
subplot(1,2,2)
dr = distanceresiduals(:,:,:);

if rewi == 1
    hold on
    errorbar(1:10,mean(dr(:,:),2),sem(dr(:,:)'),'.','color',colores(r,:))
    plot(mean(dr(:,:),2),'o','markersize',4,'markerfacecolor',colores(r,:),'markeredgecolor','none');
else
    errorbar(1:10,mean(dr(:,:),2),sem(dr(:,:)'),'.','color',colores(r,:))
    plot(mean(dr(:,:),2),'o','markersize',4,'markeredgecolor',colores(r,:),'markerfacecolor','w');
end
ylim([-5 10])
box off
axis square

xlim([0 11])
