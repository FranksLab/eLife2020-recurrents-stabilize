clearvars
close all
clc

%% Pick out files with 'kwik' in its name and put each in one cell
Catalog = 'ExperimentCatalog_bulb_awk_kx_share.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include));
Kindex = find(T.include);

%%
clear *_time
for k = 1:length(KWIKfiles)
    [efd] = loadEFD(KWIKfiles{k});
    FVs = efd.ValveTimes.FVSwitchTimesOn;
    
    inj_series = floor(T.IT(k));
    kx_start = T.FTk(k);
    kx_end = T.LTk(k);
    %%
    clear tis tkx tke
    for V = 1:size(FVs,2)
        if ~isempty(FVs{V})
            tis(V) = FVs{V}(inj_series);
            tkx(V) = FVs{V}(kx_start);
            tke(V) = FVs{V}(kx_end);
        else
            tis(V) = nan;
            tkx(V) = nan;
            tke(V) = nan;
        end
    end
    %%
    inj_time(k) = min(tis) + range(tis)*mod(T.IT(k),1);
    kxs_time(k) = min(tkx);
    kxe_time(k) = max(tke);
end

%%
onsets = (kxs_time-inj_time)/60;
offsets = (kxe_time-inj_time)/60;

%%
figure(1)
clf
printpos([100 100 200 80])
plot(mean(onsets),1,'k.','markersize',20); hold on
errbar(mean(onsets),1,sem(onsets),'k','horizontal');
plot(mean(offsets),1,'k.','markersize',20); 
errbar(mean(offsets),1,sem(offsets),'k','horizontal');
plot([0 0],[.5 1.5],'b')
set(gca,'YTick',[],'XTick',[])
box off;
xlim([-30 60])