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
    
    BinSizes = .48;
    PST = [0 .48];
    
    statelist = ['a','k'];
    
    efds = EFDloader(Catalog);
    
    [X{region},Y{region},n_bins] = CrossPopulator(Catalog, efds, BinSizes, PST, statelist); % if abs(BinSize-diff(PST))<BinSize there will be 1 bin.
    
    %%
    Xs{region} = SubsetPermuter_Cross(X{region}, Y{region}, sum(TypeStack{1}), TypeStack(1), BinSizes, n_bins);
    
    XXs{region} = cat(1,Xs{region}{:});
    YYs{region} = repmat(Y{region}{1},2,1);
end
%% no bs part.
n_trials = 7;
n_stim = 6;


for region = 1:2
    for state = 1:3
        if state <3
            sample = Xs{region}{:,:,state}(Y{region}{state}>1,:)';
            COR_real{state,region} = corr(sample,'type','Pearson');
            COR_real{state,region} = COR_real{state,region} + diag(diag(nan(length(COR_real{state,region}))));
        end
        
        if state == 3
            sample = XXs{region}(YYs{region}>1,:)';
            COR_real{state,region} = corr(sample,'type','Pearson');
            COR_real{state,region} = COR_real{state,region}(1:end/2,end/2+1:end);
        end
        %%
        clear within across
        if state <3
            for m = 1:(n_trials*n_stim)
                for n = m:(n_trials*n_stim)
                    if Y{region}{1}(m+n_trials) == Y{region}{1}(n+n_trials)
                        within{state,region}(m,n) = COR_real{state,region}(m,n);
                    else
                        across{state,region}(m,n) = COR_real{state,region}(m,n);
                    end
                end
            end
        end
        
        if state == 3
            for m = 1:(n_trials*n_stim)
                for n = 1:(n_trials*n_stim)
                    if Y{region}{1}(m+n_trials) == Y{region}{1}(n+n_trials)
                        within{state,region}(m,n) = COR_real{state,region}(m,n);
                    else
                        across{state,region}(m,n) = COR_real{state,region}(m,n);
                    end
                end
            end
        end
        w_real{state,region} = within{state,region}(within{state,region}~=0);
        a_real{state,region} = across{state,region}(across{state,region}~=0);
        
        w{state,region} = within{state,region};
        w{state,region}(w{state,region} == 0) = nan;
        
        a{state,region} = across{state,region};
        a{state,region}(a{state,region} == 0) = nan;
        
        mw_real(state,region) = nanmean(w_real{state,region});
        ma_real(state,region) = nanmean(a_real{state,region});
        wadiff_real(state,region) =  mw_real(state,region) -  ma_real(state,region);
        
        
    end
end

%% bs part.
n_trials = 7;
n_stim = 6;


for region = 1:2
    for bs = 1:1000
        samplesize = size(XXs{region},2);
                bsidx = randi(size(XXs{region},2),samplesize,1);

        for state = 1:3
            if state <3
                sample = Xs{region}{:,:,state}(Y{region}{state}>1,bsidx)';
                
                COR{state,region} = corr(sample,'type','Pearson');
                COR{state,region} = COR{state,region} + diag(diag(nan(length(COR{state,region}))));
                corstorage{state,region}(:,:,bs) = COR{state,region};
            end
            
            if state == 3
                sample = XXs{region}(YYs{region}>1,bsidx)';
                COR{state,region} = corr(sample,'type','Pearson');
                COR{state,region} = COR{state,region}(1:end/2,end/2+1:end);
                corstorage{state,region}(:,:,bs) = COR{state,region};
            end
            
            %%
            clear within across
            if state <3
                for m = 1:(n_trials*n_stim)
                    for n = m:(n_trials*n_stim)
                        if Y{region}{1}(m+n_trials) == Y{region}{1}(n+n_trials)
                            within{state,region}(m,n) = COR{state,region}(m,n);
                        else
                            across{state,region}(m,n) = COR{state,region}(m,n);
                        end
                    end
                end
            end
            
            if state == 3
                for m = 1:(n_trials*n_stim)
                    for n = 1:(n_trials*n_stim)
                        if Y{region}{1}(m+n_trials) == Y{region}{1}(n+n_trials)
                            within{state,region}(m,n) = COR{state,region}(m,n);
                        else
                            across{state,region}(m,n) = COR{state,region}(m,n);
                        end
                    end
                end
            end
            w{state,region} = within{state,region}(within{state,region}~=0);
            a{state,region} = across{state,region}(across{state,region}~=0);
            
            mw(state,bs,region) = nanmean(w{state,region});
            ma(state,bs,region) = nanmean(a{state,region});
            wadiff(state,bs,region) =  mw(state,bs,region) -  ma(state,bs,region);
        end
    end
end
%%
colores = [0.2 0.2 0.2; 0.051 0.447 0.7294; 101/255 44/255 144/255; .7 .7 .7];

%%
figure(3)
printpos([200 400 850 200])
clf

for region = 1:2
    for cc = 1:2
        subplotpos(4,1,cc+region*2-2,1,.1)
        imagesc(COR_real{cc,region});
        caxis(round(caxis*10)/10);
        axis square;
        set(gca,'XTick',[],'YTick',[],'clipping','off')
        colormap(pmkmp(128,'LinearL'));
        caxis([.2 .9]);
        ccc = caxis;
        text(21,44,[num2str(ccc(1)),' - ',num2str(ccc(2))],'FontSize',8,'color','k','horizontalalignment','center')
        title(statelist(cc))
    end
end

figure(4)
printpos([200 400 850 200])
clf
for region = 1:2
    
    subplotpos(4,1,region,1,.1)
    cc = 3;
    imagesc(COR_real{cc,region});
    caxis(round(caxis*10)/10);
    axis square;
    set(gca,'XTick',[],'YTick',[],'clipping','off')
    colormap(pmkmp(128,'LinearL'));
    caxis([0 .5]);
    ccc = caxis;
    text(21,44,[num2str(ccc(1)),' - ',num2str(ccc(2))],'FontSize',8,'color','k','horizontalalignment','center')
    
    pccs = get(gca,'Position');
    ht = pccs(3);
end

%% violins
figure(43)
printpos([500 500 400 200])
clf
% awake ob
subplot(1,4,1)
violinPlot(squeeze(mw(1,:,1))','histOri', 'left', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(1,:)}, 'xValues', [1], 'globalNorm',0,'histopt',1,'divFactor',1.2)
hold on
violinPlot(squeeze(ma(1,:,1))','histOri', 'right', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(1,:).^.3}, 'xValues', [2.2], 'globalNorm',0,'histopt',1,'divFactor',1.2)
alpha(.95)
box off;
set(gca,'XTick',[])
ylim([0.3 0.9])

% awake pcx
subplot(1,4,2)
violinPlot(squeeze(mw(1,:,2))','histOri', 'left', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(1,:)}, 'xValues', [1], 'globalNorm',0,'histopt',1,'divFactor',1.2)
hold on
violinPlot(squeeze(ma(1,:,2))','histOri', 'right', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(1,:).^.3}, 'xValues', [2.2], 'globalNorm',0,'histopt',1,'divFactor',1.2)
alpha(.95)
box off;
set(gca,'XTick',[])
ylim([0.3 0.9])

% kx ob
subplot(1,4,3)
violinPlot(squeeze(mw(2,:,1))','histOri', 'left', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(2,:)}, 'xValues', [1], 'globalNorm',0,'histopt',1,'divFactor',1.2)
hold on
violinPlot(squeeze(ma(2,:,1))','histOri', 'right', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(2,:).^.3}, 'xValues', [2.2], 'globalNorm',0,'histopt',1,'divFactor',1.2)
alpha(.95)
box off;
set(gca,'XTick',[])
ylim([0.3 0.9])

% kx pcx
subplot(1,4,4)
violinPlot(squeeze(mw(2,:,2))','histOri', 'left', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(2,:)}, 'xValues', [1], 'globalNorm',0,'histopt',1,'divFactor',1.2)
hold on
violinPlot(squeeze(ma(2,:,2))','histOri', 'right', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(2,:).^.3}, 'xValues', [2.2], 'globalNorm',0,'histopt',1,'divFactor',1.2)
alpha(.95)
box off;
set(gca,'XTick',[])
ylim([0.3 0.9])


%% separation violins
figure(44)
printpos([500 500 250 200])
clf
% ob
subplot(1,2,1)
violinPlot(squeeze(wadiff(1,:,1))','histOri', 'left', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(1,:)+[.05 0 0]}, 'xValues', [1], 'globalNorm',0,'histopt',1,'divFactor',1.2)
hold on
violinPlot(squeeze(wadiff(1,:,2))','histOri', 'right', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(1,:)}, 'xValues', [2.2], 'globalNorm',0,'histopt',1,'divFactor',1.2)
alpha(.95)
box off;
set(gca,'XTick',[])
ylim([0 .35])

subplot(1,2,2)
violinPlot(squeeze(wadiff(2,:,1))','histOri', 'left', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(2,:)+[.05 0 0]}, 'xValues', [1], 'globalNorm',0,'histopt',1,'divFactor',1.2)
hold on
violinPlot(squeeze(wadiff(2,:,2))','histOri', 'right', 'widthDiv', [1 1], 'showMM', 2, ...
    'color',  {colores(2,:)}, 'xValues', [2.2], 'globalNorm',0,'histopt',1,'divFactor',1.2)
alpha(.95)
box off;
set(gca,'XTick',[])
ylim([0 .35])

%% cross-state separation violins
figure(45)
printpos([500 500 400 200])
clf
subplot(1,2,1)
violinPlot(squeeze(wadiff(3,:,1))','histOri', 'left', 'widthDiv', [1 1], 'showMM', 2, ...
'color',  {colores(1,:)+[.05 0 0]}, 'xValues', [1], 'globalNorm',0,'histopt',1,'divFactor',1.4)
hold on
violinPlot(squeeze(wadiff(3,:,2))','histOri', 'right', 'widthDiv', [1 1], 'showMM', 2, ...
'color',  {colores(1,:)}, 'xValues', [2.2], 'globalNorm',0,'histopt',1,'divFactor',1.4)
alpha(.95)
box off;
set(gca,'XTick',[],'Clipping','off')
ylim([0 .15])

%% hypothesis testing for OB PCx differences in correlation separation.
% different-odor bulb is different from different-odor pcx
% awake diff vs awake diff, kx diff vs kx diff
% also
% null hypothesis: it doesn't matter whether pseudopopulation is built from
% ob or pcx units.

% we need Y, Xs, XXs, YYs. Y and YYs don't change so just leave it.
% we need Xs, and XXs, bootstrap scrambled versions.
% and eventually we need mw (not really), ma, and wadiff null distributions

clear w a within across

for state = 1:2
    Xmix{state} = [Xs{1}{:,:,state},Xs{2}{:,:,state}];
end

for bs = 1:1000
    for region = 1:2
        for state = 1:2
            bsidxmix{region} = randi(size(Xmix{state},2),size(XXs{region},2),1);
            Xs_bs{region}{:,:,state} = Xmix{state}(:,bsidxmix{region});
        end
        XXs_bs{region} = cat(1,Xs_bs{region}{:});
    end
    for region = 1:2
        for state = 1:3
            
            if state <3
                sample = Xs_bs{region}{:,:,state}(Y{region}{1}>1,:)';
                COR{state,region} = corr(sample,'type','Pearson');
                COR{state,region} = COR{state,region} + diag(diag(nan(length(COR{state,region}))));
            end
            
            if state == 3
                COR{state,region} = corr(XXs_bs{region}(YYs{region}>1,:)','type','Pearson');
                COR{state,region} = COR{state,region}(1:end/2,end/2+1:end);
            end
             %%
            clear within across
            if state <3
                for m = 1:(n_trials*n_stim)
                    for n = m:(n_trials*n_stim)
                        if Y{region}{1}(m+n_trials) == Y{region}{1}(n+n_trials)
                            within{state,region}(m,n) = COR{state,region}(m,n);
                        else
                            across{state,region}(m,n) = COR{state,region}(m,n);
                        end
                    end
                end
            end
            
            if state == 3
                for m = 1:(n_trials*n_stim)
                    for n = 1:(n_trials*n_stim)
                        if Y{region}{1}(m+n_trials) == Y{region}{1}(n+n_trials)
                            within{state,region}(m,n) = COR{state,region}(m,n);
                        else
                            across{state,region}(m,n) = COR{state,region}(m,n);
                        end
                    end
                end
            end
            w{state,region} = within{state,region}(within{state,region}~=0);
            a{state,region} = across{state,region}(across{state,region}~=0);
            
            mw_null(state,bs,region) = nanmean(w{state,region});
            ma_null(state,bs,region) = nanmean(a{state,region});
            wadiff_null(state,bs,region) =  mw_null(state,bs,region) -  ma_null(state,bs,region);
        end
    end
end

%%
% wadiff is higher in pcx than in bulb
for state = 1:3
    pbdiffwadiff(state) = sum((wadiff_null(state,:,2) - wadiff_null(state,:,1)) >= (wadiff_real(state,2)-wadiff_real(state,1))) / numel(wadiff_null(state,:,1));
end