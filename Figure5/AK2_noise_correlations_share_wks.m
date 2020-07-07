clearvars
close all
clc

specialparams.FRlim = 1/100;
specialparams.UVlim = 50;
specialparams.DFRlim = 100;

statelist = ['a','k'];

Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_pcx_awk_kx_F.txt';
[TypeIdx_p, ~] = CellTyper (Catalog, 'Stable', specialparams);

Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_bulb_awk_kx_F_sim.txt';
[TypeIdx_b, ~] = CellTyper (Catalog, 'Stable', specialparams);


binsz = .03;
stepsz = .015;
winsz = .03;
winsteps = [-.2:stepsz:.4];


for bs = 1:length(winsteps)
    PST = [winsteps(bs) winsteps(bs)+winsz+.001];
    
    %% pcx trial stack population vectors
    Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_pcx_awk_kx_F.txt';
    
    T = readtable(Catalog, 'Delimiter', ' ');
    KWIKfiles = T.kwikfile(logical(T.include));
    Kindex = find(T.include);
    
    efds = EFDloader(Catalog);
    [Xp,Y,n_bins] = CrossPopulator_single_experiments(Catalog, efds, binsz, PST, statelist); % if abs(BinSize-diff(PST))<BinSize there will be 1 bin.
    
    %% bulb pcx trial stack population vectors
    Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_bulb_awk_kx_F_sim.txt';
    
    T = readtable(Catalog, 'Delimiter', ' ');
    KWIKfiles = T.kwikfile(logical(T.include));
    Kindex = find(T.include);
    
    statelist = ['a','k'];
    
    efds = EFDloader(Catalog);
    [Xb,Y,n_bins] = CrossPopulator_single_experiments(Catalog, efds, binsz, PST, statelist); % if abs(BinSize-diff(PST))<BinSize there will be 1 bin.
    
    %%
    clear NC
    %%
    for k = 1:length(efds)
        Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_pcx_awk_kx_F.txt';
        T = readtable(Catalog, 'Delimiter', ' ');
        KWIKfiles = T.kwikfile(logical(T.include));
        [SpikeTimes] = SpikeTimesPro(FindFilesKK(KWIKfiles{k}));
        [ypos{k},pttime,asym] = WaveStats(SpikeTimes.Wave);
        % add back an empty spot to use celltyper
        ypos{k} = [nan; ypos{k}]; ypos{k} = ypos{k}(TypeIdx_p{k,1});
        %%
        for state = 1:2
            for stim = 2:7
                stacko = [Xb{k,state}(Y==stim,TypeIdx_b{k}), Xp{k,state}(Y==stim,TypeIdx_p{k})];
                b_cells = 1:size(Xb{k,state}(:,TypeIdx_b{k}),2);
                p_cells = (size(Xb{k,state}(:,TypeIdx_b{k}),2)+1):(size(Xb{k,state}(:,TypeIdx_b{k}),2)+size(Xp{k,state}(:,TypeIdx_p{k}),2));
                
                % z scoring the stacko because of miura
                s = bsxfun(@minus,stacko,mean(stacko));
                stacko = bsxfun(@rdivide,s,std(stacko));
                
                c = corr(stacko);
                c = c + diag(diag(nan(length(c))));
                NC(k,state,stim-1,1) = nanmean(reshape(c(b_cells,b_cells),[],1));
                NC(k,state,stim-1,2) = nanmean(reshape(c(p_cells,p_cells),[],1));
                NC(k,state,stim-1,3) = nanmean(reshape(c(b_cells,p_cells),[],1));
                
                NCdist{bs}{k,state,stim-1,1} = (reshape(c(b_cells,b_cells),[],1));
                NCdist{bs}{k,state,stim-1,2} = (reshape(c(p_cells,p_cells),[],1));
                NCdist{bs}{k,state,stim-1,3} = (reshape(c(b_cells,p_cells),[],1));
                
                BPmn{bs,k,state}(stim-1,:) = nanmean(c(b_cells,p_cells));
                BPmx{bs,k,state}(stim-1,:) = max(c(b_cells,p_cells));
                
                PPmn{bs,k,state}(stim-1,:) = nanmean(c(p_cells,p_cells));
                PPmx{bs,k,state}(stim-1,:) = max(c(p_cells,p_cells));
            end
        end
    end
    
    NCall{bs} = squeeze(mean(NC,3));
    
end

%% rearrange NCdist
X = cat(5,NCdist{:});
clear Ndr
for state = 1:2
    for conn_type = 1:3
        for bs = 1:length(winsteps)
            Ndr{state,conn_type}(bs,:) = reshape(cell2mat(X(:,state,:,conn_type,bs)),[],1);
        end
    end
end


%%
figure(3)
clf
printpos([100 100 800 300])
NCmat = cat(4,NCall{:});
titles = {'b-b','p-p','b-p'};

for conn_type = 1:3
    subplot(1,3,conn_type)
    boundedline(winsteps+winsz/2,squeeze(nanmean(NCmat(:,1,conn_type,:))),squeeze(nansem(NCmat(:,1,conn_type,:))),'k')

    hold on
    boundedline(winsteps+winsz/2,squeeze(nanmean(NCmat(:,2,conn_type,:))),squeeze(nansem(NCmat(:,2,conn_type,:))))
    
    axis square; box off; set(gca,'Clipping','off')
    xlim([winsteps(1) winsteps(end)+winsz])
    title(titles{conn_type})
    
    ylim([0 0.3])
    
    if conn_type ==3
        ylim([-.02 0.06])
    end
    
end
