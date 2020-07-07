function [X,Y,n_bins,trd] = CrossPopulatorFifteenTrials_ntng(Catalog, efds, BinSizes, PST, statelist)

%% Read that catalog
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include));
Kindex = find(T.include);

for state = 1:length(statelist)
    FT = T.(['FT',statelist(state)]);
    FT = FT(Kindex);
    LT = T.(['LT',statelist(state)]);
    LT = LT(Kindex);
    
    %% Concatenate the catalogued data [X,Y,n_bins]
    
    % Preallocate
    traindata = cell(length(KWIKfiles),length(BinSizes));
    n_bins = nan(size(BinSizes));
    
    % Collect binned data, flattened into a column
    for R = 1:length(KWIKfiles)
        efd = efds(R);
        Trials = 1:15;
        VOI = 2:11;
        
        [Raster_S] = efd.ValveSpikes.RasterAlign(:,1,:);
        Raster_S = squeeze(Raster_S(:,1,:));
        
        %%
        for BS = 1:length(BinSizes)
            [Y{state},traindata{R,BS},n_bins(BS)] = BinRearranger(Raster_S([1,VOI],:,:),PST,BinSizes(BS),Trials);
            trd{R,state} = traindata{R,BS};
            
        end
    end
    
    % Concatenate across experiments.
    for BS = 1:length(BinSizes)
        X{state,BS} = cat(2,traindata{:,BS});
    end
end

end