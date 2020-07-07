clearvars
close all
clc


for region = 1:2
    %% Pick out files with 'kwik' in its name and put each in one cell
    if region == 1
        Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_bulb_awk_kx_F.txt';
    else
        Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_pcx_awk_kx_F.txt';
    end
    
    T = readtable(Catalog, 'Delimiter', ' ');
    KWIKfiles = T.kwikfile(logical(T.include));
    Kindex = find(T.include);
    
    specialparams.FRlim = 1/100;
    specialparams.UVlim = 50;
    specialparams.DFRlim = 100;
    
    [TypeIdx, TypeStack] = CellTyper (Catalog, 'Stable', specialparams);
    
    %%
    PST = [0, 0.5];
    
    for k = 1:length(KWIKfiles)
        clear KDF
        efd(k) = loadEFD(KWIKfiles{k});
        
        if strcmp(T.VOI(Kindex(k)),'A')
            VOI = [4,7,8,12,15,16];
        elseif strcmp(T.VOI(Kindex(k)),'C')
            VOI = [6,7,8,10,11,12];
        end
        
        TOI{1} = (T.FTa(Kindex(k))):T.LTa(Kindex(k));
        TOI{2} = (T.FTk(Kindex(k))):T.LTk(Kindex(k));
        for state = 1:2
            Trials = TOI{state};
            [Scores(state),~] = SCOmaker_NoBlank(KWIKfiles{k},{Trials});
        end
        CrossState{1}{k} = Scores(1).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(1).auROC(VOI,TypeIdx{k,1},1) > .5 & ~(Scores(2).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(2).auROC(VOI,TypeIdx{k,1},1) > .5);
        CrossState{2}{k} = Scores(2).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(2).auROC(VOI,TypeIdx{k,1},1) > .5 & ~(Scores(1).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(1).auROC(VOI,TypeIdx{k,1},1) > .5);
        CrossState{3}{k} = Scores(1).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(1).auROC(VOI,TypeIdx{k,1},1)>.5 & Scores(2).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(2).auROC(VOI,TypeIdx{k,1},1)>.5;
        
        FracPres{region}(k) = sum(CrossState{3}{k}(:)) / (sum(CrossState{3}{k}(:)) + sum(CrossState{1}{k}(:)));
        NumCells{region}(k) = size(CrossState{1}{k},2);

    end
    
    for c = 1:3
        CS{c,region} = reshape(cell2mat(CrossState{c}),[],1);
    end
    
    pcts_real{region}(1,1) = (sum(CS{1,region})) / ((sum(CS{1,region})) + (sum(CS{3,region})));
    pcts_real{region}(1,2) = (sum(CS{3,region})) / ((sum(CS{1,region})) + (sum(CS{3,region})));
    pcts_real{region}(2,1) = (sum(CS{2,region})) / ((sum(CS{2,region})) + (sum(CS{3,region})));
    pcts_real{region}(2,2) = (sum(CS{3,region})) / ((sum(CS{2,region})) + (sum(CS{3,region})));
    
    %% bootstrap

    for bs = 1:1000
        bsidx = randi(size(CS{1,region},1),size(CS{1,region},1),1);
        pcts_bs{region}(bs,1,1) = (sum(CS{1,region}(bsidx))) / ((sum(CS{1,region}(bsidx))) + (sum(CS{3,region}(bsidx))));
        pcts_bs{region}(bs,1,2) = (sum(CS{3,region}(bsidx))) / ((sum(CS{1,region}(bsidx))) + (sum(CS{3,region}(bsidx))));
        pcts_bs{region}(bs,2,1) = (sum(CS{2,region}(bsidx))) / ((sum(CS{2,region}(bsidx))) + (sum(CS{3,region}(bsidx))));
        pcts_bs{region}(bs,2,2) = (sum(CS{3,region}(bsidx))) / ((sum(CS{2,region}(bsidx))) + (sum(CS{3,region}(bsidx))));
    end
    
    ll{region} = squeeze(prctile(pcts_bs{region},2.5));
    ul{region} = squeeze(prctile(pcts_bs{region},97.5));

end
%% bars
figure(1)
printpos([100 400 100 220])
clf
colores = [1 .2 .2; 0.4 0.4 0.4];

% state,cc ->
% pct_real{region}(1,1) % awake-only / awake-only + robust
% pct_real{region}(1,2) % robust / awake-only + robust

subplot(2,1,1)
for region = [2,1]
    state = 1; cc = 2;
    bar(region,100*pcts_real{region}(state,cc),'facecolor',colores(region,:))
    hold on; box off; set(gca,'XTick',[],'clipping','off'); xlim([0 3]); ylim([0 100])
    errorbar(region,100*pcts_real{region}(state,cc),100*pcts_real{region}(state,cc)-100*ll{region}(state,cc),100*ul{region}(state,cc)-100*pcts_real{region}(state,cc),'k');
end
ylim([0 50])
axis square


%% box and whisker
figure(2)
printpos([100 400 100 220])
clf
colores = [1 .2 .2; 0.4 0.4 0.4; 1 .2 .2; 0.4 0.4 0.4];

% state,cc ->
% pct_real{region}(1,1) % awake-only / awake-only + robust
% pct_real{region}(1,2) % robust / awake-only + robust
% pct_real{region}(2,2) % robust / kx-only + robust


subplot(1,1,1)
state = 1; cc = 2;
figure_boxplot([100*pcts_bs{1}(:,1,2), 100*pcts_bs{2}(:,1,2),100*pcts_bs{1}(:,2,2), 100*pcts_bs{2}(:,2,2)],{[];[];[];[]},[],[],'horizontal',colores(:,:));

hold on; box off; set(gca,'XTick',[],'clipping','off'); xlim([0 5]); ylim([0 70])


%% actual mean difference

% null hypothesis for significance test
CS1 = cat(1,CS{1,:});
CS2 = cat(1,CS{2,:});
CS3 = cat(1,CS{3,:});

clear mdnh
for bs = 1:1000
    for region= 1:2
        bsidx = randi(size(CS1,1),size(CS{1,region},1),1);
        CS1bs = CS1(bsidx);
        CS1bs = CS2(bsidx);
        CS3bs = CS3(bsidx);
        
        pct_preserved(region) = sum(CS3bs) / (sum(CS1bs) + sum(CS3bs)); % for testing awk->kx preservation
%         pct_preserved(region) = sum(CS3bs) / (sum(CS1bs) + sum(CS3bs)); % for testing kx->awk preservation
    end
    mdnh(bs) = pct_preserved(2) - pct_preserved(1);
end
   
md = pcts_real{2}(1,2) - pcts_real{1}(1,2); % for testing awk->kx preservation
% md = pcts_real{2}(2,2) - pcts_real{1}(2,2); % for testing kx->awk preservation


pp = sum(mdnh > md)/numel(mdnh);
