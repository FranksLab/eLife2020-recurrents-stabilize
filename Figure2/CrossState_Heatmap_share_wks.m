clearvars
close all
clc

%% choose a set of recordings
% Catalog = 'ExperimentCatalog_bulb_awk_kx_share.txt';
Catalog = 'ExperimentCatalog_pcx_awk_kx_share.txt';

specialparams.FRlim = 1/100;
specialparams.UVlim = 50;
specialparams.DFRlim = 100;

[TypeIdx, TypeStack] = CellTyper (Catalog, 'Stable', specialparams);

T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include));
Kindex = find(T.include);

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
    
    if size(Scores(1).auROC,2) ~= size(TypeIdx{k,1},1)
        disp('CellTyper error. Different numbers of units indexed')
        break
    end
    CrossState{1}{k} = Scores(1).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(1).auROC(VOI,TypeIdx{k,1},1) > .5 & ~(Scores(2).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(2).auROC(VOI,TypeIdx{k,1},1) > .5);
    CrossState{2}{k} = Scores(2).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(2).auROC(VOI,TypeIdx{k,1},1) > .5 & ~(Scores(1).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(1).auROC(VOI,TypeIdx{k,1},1) > .5);
    CrossState{3}{k} = Scores(1).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(1).auROC(VOI,TypeIdx{k,1},1)>.5 & Scores(2).AURp(VOI,TypeIdx{k,1},1) < .05 & Scores(2).auROC(VOI,TypeIdx{k,1},1)>.5;
    
    %%
    [SpikeTimes] = SpikeTimesPro(FindFilesKK(KWIKfiles{k}));
    [ypos,pttime,asym] = WaveStats(SpikeTimes.Wave);
    ypos = [nan; ypos]; ypos = ypos(TypeIdx{k,1});
    
    % more negative is more superficial. lower index is more superficial.
    % we want ventral on top. ventral is more superficial. we want lower index
    % on top.
    [sortpos,posidx] = sort(ypos);
    
    %%
    figure(k)
    clf
    printpos([100 100 250 sum(TypeIdx{k,1})*6])
    colores = [1,1,1; 0.75 0.75 0.75; 0.491 0.797 0.88; 180/255 34/255 180/255]; % white scheme magenta
    
    % white scheme
    CT = flipud(cbrewer('div','RdBu',64));
    CT = CT(3:end-3,:);
    
    subplotpos(3,1,3,1,.3)
    CSmap = zeros(size(CrossState{1}{k}));
    for cs = 1:3
        CSmap = CSmap + cs*CrossState{cs}{k};
    end
    ax3 = gca;
    imagesc(CSmap(:,posidx)'); caxis([0 4])
    set(gca,'YTick',get(gca,'YLim')+[.5 -.5],'XTick',[])
    
    colormap(ax3,colores);
    text(7,-1,num2str(k))
    
    for state = 1:2
        subplotpos(3,1,state,1,.3)
        RImap = Scores(state).auROC(VOI,TypeIdx{k,1},1);
        
        ax{state} = gca;
        imagesc(RImap(:,posidx)'); caxis([0 1])
        set(gca,'YTick',get(gca,'YLim')+[.5 -.5],'XTick',[])
        
        colormap(ax{state},CT);
    end
    
end

