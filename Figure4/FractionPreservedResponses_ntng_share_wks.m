clearvars
close all
clc

%% Pick out files with 'kwik' in its name and put each in one cell
Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_KX_Ntng_F.txt';

T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include));
Kindex = find(T.include);

specialparams.FRlim = 1/100;
specialparams.UVlim = 50;
specialparams.DFRlim = 100;

[TypeIdx, TypeStack] = CellTyper (Catalog, 'StableNTNG', specialparams);
%%
for k = 1:length(KWIKfiles)
    clear KDF
    efd(k) = loadEFD(KWIKfiles{k});
    VOI = [2 3 4 5 7 11];
    TOI{1} = ((T.FTa(Kindex(k))):(T.LTa(Kindex(k))));
    TOI{2} = (T.FTk(Kindex(k))):(T.LTk(Kindex(k)));
    
    for state = 1:2
        Trials = TOI{state};
        [Scores(state),~] = SCOmaker_Beast_NoBlank(KWIKfiles{k},{Trials});
    end
    for celltype = 1:2
        
        alpha = .05;
        CrossState{1,celltype}{k} = squeeze(Scores(1).AURp(VOI,1,TypeIdx{k,celltype},2) < alpha & Scores(1).auROC(VOI,1,TypeIdx{k,celltype},2) > .5 & ~(Scores(2).AURp(VOI,1,TypeIdx{k,celltype},2) < alpha & Scores(2).auROC(VOI,1,TypeIdx{k,celltype},2) > .5));
        CrossState{2,celltype}{k} = squeeze(Scores(2).AURp(VOI,1,TypeIdx{k,celltype},2) < alpha & Scores(2).auROC(VOI,1,TypeIdx{k,celltype},2) > .5 & ~(Scores(1).AURp(VOI,1,TypeIdx{k,celltype},2) < alpha & Scores(1).auROC(VOI,1,TypeIdx{k,celltype},2) > .5));
        CrossState{3,celltype}{k} = squeeze(Scores(1).AURp(VOI,1,TypeIdx{k,celltype},2) < alpha & Scores(1).auROC(VOI,1,TypeIdx{k,celltype},2)>.5 & Scores(2).AURp(VOI,1,TypeIdx{k,celltype},2) < alpha & Scores(2).auROC(VOI,1,TypeIdx{k,celltype},2)>.5);
    end
end

%%
for celltype = 1:2
    
    %% bootstrap over responses
    for c = 1:3
        CS{c,celltype} = reshape(cell2mat(CrossState{c,celltype}),[],1);
    end
    
    pcts_real{celltype}(1,1) = (sum(CS{1,celltype})) / ((sum(CS{1,celltype})) + (sum(CS{3,celltype})));
    pcts_real{celltype}(1,2) = (sum(CS{3,celltype})) / ((sum(CS{1,celltype})) + (sum(CS{3,celltype})));
    pcts_real{celltype}(2,1) = (sum(CS{2,celltype})) / ((sum(CS{2,celltype})) + (sum(CS{3,celltype})));
    pcts_real{celltype}(2,2) = (sum(CS{3,celltype})) / ((sum(CS{2,celltype})) + (sum(CS{3,celltype})));
    %
    
    for bs = 1:10%00
        bsidx = randi(size(CS{1,celltype},1),size(CS{1,celltype},1),1);
        pcts_bs{celltype}(bs,1,1) = (sum(CS{1,celltype}(bsidx))) / ((sum(CS{1,celltype}(bsidx))) + (sum(CS{3,celltype}(bsidx))));
        pcts_bs{celltype}(bs,1,2) = (sum(CS{3,celltype}(bsidx))) / ((sum(CS{1,celltype}(bsidx))) + (sum(CS{3,celltype}(bsidx))));
        pcts_bs{celltype}(bs,2,1) = (sum(CS{2,celltype}(bsidx))) / ((sum(CS{2,celltype}(bsidx))) + (sum(CS{3,celltype}(bsidx))));
        pcts_bs{celltype}(bs,2,2) = (sum(CS{3,celltype}(bsidx))) / ((sum(CS{2,celltype}(bsidx))) + (sum(CS{3,celltype}(bsidx))));
    end
    
    ll{celltype} = squeeze(prctile(pcts_bs{celltype},2.5));
    ul{celltype} = squeeze(prctile(pcts_bs{celltype},97.5));
end

%% box and whisker plots
figure(2)
printpos([100 400 100 220])
clf
colores = [0.2 0.2 0.2;101/255 44/255 144/255; 0.051 0.447 0.7294; 101/255 44/255 144/255; .7 .7 .7];
colores = [1 0 0; .2 .2 .2; .7 .7 .7];
% state,cc ->
% pct_real{celltype}(1,1) % awake-only / awake-only + robust
% pct_real{celltype}(1,2) % robust / awake-only + robust

subplot(2,1,1)
state = 1; cc = 2;
figure_boxplot([100*pcts_bs{1}(:,state,cc), 100*pcts_bs{2}(:,state,cc)],{[];[]},[],[],'horizontal',colores(:,:));

hold on; box off; set(gca,'XTick',[],'clipping','off'); xlim([0 3]); ylim([0 100])
ylim([0 60])
axis square

%% actual mean difference
md = pcts_real{2}(1,2) - pcts_real{1}(1,2);

% null hypothesis for significance test
CS1 = cat(1,CS{1,:});
CS3 = cat(1,CS{3,:});
clear mdnh
for bs = 1:1000
    for region= 1:2
        bsidx = randi(size(CS1,1),size(CS{1,region},1),1);
        CS1bs = CS1(bsidx);
        CS3bs = CS3(bsidx);
        
        pct_preserved(region) = sum(CS3bs) / (sum(CS1bs) + sum(CS3bs));
    end
    mdnh(bs) = pct_preserved(2) - pct_preserved(1);
end

pp = sum(mdnh > md)/numel(mdnh);

%% heat maps
k = 14;
celltype = 2;

for celltype = 1:2
    
    figure(celltype)
    clf
    printpos([100 100 250 sum(TypeIdx{k,celltype})*6])
    colores = [1,1,1; 0.75 0.75 0.75; 0.491 0.797 0.88; 180/255 34/255 180/255]; % white scheme magenta lighten state-specifics
    
    % white scheme
    CT = flipud(cbrewer('div','RdBu',64));
    CT = CT(3:end-3,:);
    
    
    subplotpos(3,1,3,1,.3)
    CSmap = zeros(size(CrossState{1,celltype}{k}));
    for cs = 1:3
        CSmap = CSmap + cs*CrossState{cs,celltype}{k};
    end
    ax3 = gca;
    imagesc(CSmap(:,:)'); caxis([0 4])
    set(gca,'YTick',get(gca,'YLim')+[.5 -.5],'XTick',[])
    
    colormap(ax3,colores);
    text(7,-1,num2str(k))
    
    for state = 1:2
        subplotpos(3,1,state,1,.3)
        RImap = Scores(state).auROC(VOI,1,TypeIdx{k,celltype},2);
        
        ax{state} = gca;
        imagesc(RImap(:,:)'); caxis([0 1])
        set(gca,'YTick',get(gca,'YLim')+[.5 -.5],'XTick',[])
        
        colormap(ax{state},CT);
    end
    
end