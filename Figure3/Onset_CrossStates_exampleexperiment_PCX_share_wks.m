clearvars
close all
clc

%% Pick out files with 'kwik' in its name and put each in one cell
Catalog = 'ExperimentCatalog_pcx_awk_kx_share.txt';

T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include));
Kindex = find(T.include);

specialparams.FRlim = 1/100;
specialparams.UVlim = 50;
specialparams.DFRlim = 100;

[TypeIdx, TypeStack] = CellTyper (Catalog, 'Stable', specialparams);

%%
PST = [-2, 0.55];

k = 6;

[SpikeTimes] = SpikeTimesPro(FindFilesKK(KWIKfiles{k}));
[ypos,pttime,asym] = WaveStats(SpikeTimes.Wave);

% add back an empty spot to use celltyper
ypos = [nan; ypos]; ypos = ypos(TypeIdx{k,1});
pttime = [nan; pttime]; pttime = pttime(TypeIdx{k,1});
asym = [nan; asym]; asym = asym(TypeIdx{k,1});

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

for tset = 1:2
    % for finding peaks
    BinSize = 0.010;
    [KDF{tset}, ~, KDFt, KDFe] = KDFmaker(efd(k).ValveSpikes.RasterAlign(VOI,TypeIdx{k,1}), PST, BinSize, TOI{tset});
    realPST = KDFt>=PST(1) & KDFt<=PST(2);
    KDFt = KDFt(realPST);
    
    for valve = 1:(size(KDF{tset},1))
        for unit = 1:size(KDF{tset},2)
            if ~isempty(KDF{tset}{valve,unit}) && ~isempty(KDF{tset}{valve,unit})
                KDF{tset}{valve,unit} = KDF{tset}{valve,unit}(realPST);
                [mx{tset}{k}(valve,unit),mxt_temp] = max(KDF{tset}{valve,unit});
                mxt{tset}{k}(valve,unit) = KDFt(mxt_temp);
                mxd{tset}{k}(valve,unit) = FWHM(KDFt,KDF{tset}{valve,unit},.5);
                depth{tset}{k}(valve,unit) = ypos(unit);
                pt{tset}{k}(valve,unit) = pttime(unit);
                as{tset}{k}(valve,unit) = asym(unit);
                
                baseline = KDF{tset}{valve,unit}(KDFt<0);
                first_below = find(KDF{tset}{valve,unit}(KDFt>=0) <= (mean(baseline) + 2.5*std(baseline)),1);
                if ~isempty(first_below)
                    firstbelowtime = KDFt(find(KDFt>=0,1) + first_below -1);
                    onset_idx = find(KDF{tset}{valve,unit}(KDFt>=firstbelowtime) > (mean(baseline) + 2.5*std(baseline)),1);
                    if ~isempty(onset_idx)
                        onset{tset}{k}(valve,unit) = KDFt(find(KDFt>=firstbelowtime,1) + onset_idx -1);
                    else
                        onset{tset}{k}(valve,unit) = nan;
                    end
                else
                    onset{tset}{k}(valve,unit) = nan;
                end
                
            else
                KDF{tset}{valve,unit}  = nan(1,length(realPST));
                mx{tset}{k}(valve,unit) = nan;
                mxt{tset}{k}(valve,unit) = nan;
                mxd{tset}{k}(valve,unit) = nan;
                depth{tset}{k}(valve,unit) = nan;
                pt{tset}{k}(valve,unit) = nan;
                as{tset}{k}(valve,unit) = nan;
                onset{tset}{k}(valve,unit) = nan;
            end
        end
    end
end


%%
Ka = KDF{1}; Ka = Ka(:);
Kk = KDF{2}; Kk = Kk(:);

La = cat(2,onset{1}{:}); La = La(:);
Lk = cat(2,onset{2}{:}); Lk = Lk(:);
Nopks_a = La < (PST(1)+.001) | La > (PST(2)-.001);
Nopks_k = Lk < (PST(1)+.001) | Lk > (PST(2)-.001);
La(Nopks_a) = nan;
Lk(Nopks_k) = nan;

Pa = cat(2,mx{1}{:}); Pa = Pa(:);
Pk = cat(2,mx{2}{:}); Pk = Pk(:);
Pa(Nopks_a) = nan;
Pk(Nopks_k) = nan;

Da = cat(2,mxd{1}{:}); Da = Da(:);
Dk = cat(2,mxd{2}{:}); Dk = Dk(:);
Da(Nopks_a) = nan;
Dk(Nopks_k) = nan;

Ya = cat(2,depth{1}{:}); Ya = Ya(:);
Yk = cat(2,depth{2}{:}); Yk = Yk(:);
Ya(Nopks_a) = nan;
Yk(Nopks_k) = nan;

Ta = cat(2,pt{1}{:}); Ta = Ta(:);
Tk = cat(2,pt{2}{:}); Tk = Tk(:);
Ta(Nopks_a) = nan;
Tk(Nopks_k) = nan;

Aa = cat(2,as{1}{:}); Aa = Aa(:);
Ak = cat(2,as{2}{:}); Ak = Ak(:);
Aa(Nopks_a) = nan;
Ak(Nopks_k) = nan;


CSa = cat(2,CrossState{1}{:}); CSa = logical(CSa(:));
CSk = cat(2,CrossState{2}{:}); CSk = logical(CSk(:));
CSc = cat(2,CrossState{3}{:}); CSc = logical(CSc(:));

%%
figure(5)
clf
printpos([300 100 150 600])

La = cat(2,onset{1}{:}); La = La(:);
Lk = cat(2,onset{2}{:}); Lk = Lk(:);

subplot(1,1,1)
[Y,I] = sort(La,'descend');

I = I(~isnan(Y));

spacing = 50;
xleft = -.05;

for R = length(I):-1:1
    hold on
    if CSa(I(R))
        plot(KDFt(KDFt>=xleft),R*spacing+(Ka{I(R)}(KDFt>=xleft) - mean(Ka{I(R)}(KDFt<0))),'k')
        plot(La(I(R)),R*spacing+(Ka{I(R)}(find(KDFt>La(I(R)),1)) - mean(Ka{I(R)}(KDFt<0))),'.','color','k','markersize',8)
    elseif CSc(I(R))
        plot(KDFt(KDFt>=xleft),R*spacing+(Ka{I(R)}(KDFt>=xleft) - mean(Ka{I(R)}(KDFt<0))),'color','m')
        plot(La(I(R)),R*spacing+(Ka{I(R)}(find(KDFt>La(I(R)),1)) - mean(Ka{I(R)}(KDFt<0))),'.','color','m','markersize',8)
    else
        plot(KDFt(KDFt>=xleft),R*spacing+(Ka{I(R)}(KDFt>=xleft) - mean(Ka{I(R)}(KDFt<0))),'color',[.8 .8 .8])
        if ~isempty(find(KDFt>La(I(R)),1))
            plot(La(I(R)),R*spacing+(Ka{I(R)}(find(KDFt>La(I(R)),1)) - mean(Ka{I(R)}(KDFt<0))),'.','color',[.8 .8 .8],'markersize',8)
        end
    end
    hold on
    box off
    xlim([xleft .55]); ylim ([0 (length(I)+2)*spacing])
end
set(gca,'YTick',[],'XTick',[-.5 0 .5],'clipping','off')

hold on
plot([0 0],[0 100],'k')