clearvars
close all
clc

%% For plotting 'tuning' measures, percent responsive, sparsensess, sorted tuning, different cell layers
% turn on Essentials

%% Pick out files with 'kwik' in its name and put each in one cell
% Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_bulb_awk_kx_F.txt';
Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_pcx_awk_kx_F.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include));
Kindex = find(T.include);

specialparams.FRlim = 1/100;
specialparams.UVlim = 50;
specialparams.DFRlim = 100;

[TypeIdx, TypeStack] = CellTyper (Catalog, 'Stable', specialparams);
%%

% clear Tuning*
for k = 1:length(KWIKfiles);
    TOI{1} = T.FTa(Kindex(k)):T.LTa(Kindex(k));
    TOI{2} = T.FTk(Kindex(k)):T.LTk(Kindex(k));
    
    if strcmp(T.VOI(Kindex(k)),'A')
        VOI = [4,7,8,12,15,16];
    elseif strcmp(T.VOI(Kindex(k)),'C')
        VOI = [6,7,8,10,11,12];
    end
    
    for state = 1:2
        Trials = TOI{state};
        [Scores,~] = SCOmaker_NoBlank(KWIKfiles{k},{Trials});
        
        
        pr{k,state} = Scores.auROC(VOI,TypeIdx{k,1},1)>.5 & Scores.AURp(VOI,TypeIdx{k,1},1)<.05;
        nr{k,state} = Scores.auROC(VOI,TypeIdx{k,1},1)<.5 & Scores.AURp(VOI,TypeIdx{k,1},1)<.05;
    end
end
prs{1} = cat(2,pr{:,1});
prs{2} = cat(2,pr{:,2});
nrs{1} = cat(2,nr{:,1});
nrs{2} = cat(2,nr{:,2});


%% prs and nrs bootstrap

% to get confidence intervals
for bs = 1:1000
    bsidx = randi(size(prs{1},2),size(prs{1},2),1);
    for state = 1:2
        prs_pct(state,bs) = 100 * sum(sum(prs{state}(:,bsidx))) / numel(prs{state}(:,bsidx));
        nrs_pct(state,bs) = 100 * sum(sum(nrs{state}(:,bsidx))) / numel(nrs{state}(:,bsidx));
    end
end

% actual mean difference
prs_md = 100 * abs((sum(sum(prs{2})) / numel(prs{2})) - (sum(sum(prs{1})) / numel(prs{1})));
nrs_md = 100 * abs((sum(sum(nrs{2})) / numel(nrs{2})) - (sum(sum(nrs{1})) / numel(nrs{1})));

prs_mix = cell2mat(prs);
nrs_mix = cell2mat(nrs);

% null hypothesis for significance test - what if awake and kx cells were
% in the same distribution?
for bs = 1:1000
    bsidx_a = randi(size(prs_mix,2),size(prs{1},2),1);
    bsidx_b = randi(size(prs_mix,2),size(prs{1},2),1);
    
    prs_nh(bs) = 100 * abs((sum(sum(prs_mix(:,bsidx_a))) / numel(prs_mix(:,bsidx_a))) - (sum(sum(prs_mix(:,bsidx_b))) / numel(prs_mix(:,bsidx_b))));
    nrs_nh(bs) = 100 * abs((sum(sum(nrs_mix(:,bsidx_a))) / numel(nrs_mix(:,bsidx_a))) - (sum(sum(nrs_mix(:,bsidx_b))) / numel(nrs_mix(:,bsidx_b))));
end

prs_pp = sum(prs_nh > prs_md)/numel(prs_nh);
nrs_pp = sum(nrs_nh > nrs_md)/numel(nrs_nh);

%% Raw rates for sparseness calculation
for k = 1:length(KWIKfiles);
    TOI{1} = T.FTa(Kindex(k)):T.LTa(Kindex(k));
    TOI{2} = T.FTk(Kindex(k)):T.LTk(Kindex(k));
    for state = 1:2
        Trials = TOI{state};
        [Scores,~] = SCOmaker_NoBlank(KWIKfiles{k},{Trials});
        RR{k,state} = Scores.RawRate(VOI,TypeIdx{k,1},1);
    end
end

%% sparseness bootstrap

% to get confidence intervals
for bs = 1:1000
    bsidx = randi(size(cell2mat(RR(:,1)'),2),size(cell2mat(RR(:,1)'),2),1);
    for state = 1:2
        RRs = cell2mat(RR(:,state)');
        SparseVar = RRs(:,bsidx);
        [L_bs, P_bs] = Sparseness(SparseVar);
        SL_bs(state,bs) = nanmean(L_bs);
        SP_bs(state,bs) = nanmean(P_bs);
    end
end

% actual mean difference
for state = 1:2
    RRs = cell2mat(RR(:,state)');
    SparseVar = RRs(:,:);
    [L(state,:),P(state,:)] = Sparseness(SparseVar);
end
L_md = abs(diff(nanmean(L')));
P_md = abs(diff(nanmean(P')));

% null hypothesis for significance test
RRmix = cat(2,RR{:});
for bs = 1:1000
    for state = 1:2
        bsidx = randi(size(RRmix,2),size(cell2mat(RR(:,1)'),2),1);
        RRs = RRmix(:,bsidx);
        SparseVar = RRs(:,:);
        [L_nh(state,:),P_nh(state,:)] = Sparseness(SparseVar);
    end
    L_mdnh(bs) = abs(diff(nanmean(L_nh')));
    P_mdnh(bs) = abs(diff(nanmean(P_nh')));
end

L_pp = sum(L_mdnh > L_md)/numel(L_mdnh);
P_pp = sum(P_mdnh > P_md)/numel(P_mdnh);

%% box and whisker plots
colores = [0.2 0.2 0.2; 0.051 0.447 0.7294];

figure(3)
clf
printpos([200 200 500 150])
subplot(1,4,1)
figure_boxplot(prs_pct(:,:)',{[];[]},[],[],'horizontal',colores);

box off; ylim([0 18]); xlim([0 3]);
set(gca,'XTick',[1.5],'XTickLabel',[],'Clipping','off');
set(gca,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1.2)
ylabel('cell-odor pairs (%)')

subplot(1,4,2)
figure_boxplot(nrs_pct(:,:)',{[];[]},[],[],'horizontal',colores);

box off; ylim([0 18]); xlim([0 3]);
set(gca,'XTick',[1.5],'XTickLabel',[],'Clipping','off');
set(gca,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1.2)
ylabel('cell-odor pairs (%)')

subplot(1,4,3)
figure_boxplot(SL_bs(:,:)',{[];[]},[],[],'horizontal',colores);
box off; ylim([0 1]);  xlim([0 3]);
set(gca,'XTick',[1.5],'XTickLabel',{'lifetime'},'YTick',[0 1],'Clipping','off');
set(gca,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1.2)
ylabel('sparseness')

subplot(1,4,4)
figure_boxplot(SP_bs(:,:)',{[];[]},[],[],'horizontal',colores);
box off; ylim([0 1]);  xlim([0 3]);
set(gca,'XTick',[1.5],'XTickLabel',{'pop.'},'YTick',[0 1],'Clipping','off');
set(gca,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1.2)
ylabel('sparseness')



