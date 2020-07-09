function [X,Y,n_bins] = PseudoPopulator_loadingSets_ntng(Catalog, efds, BinSizes, PST)

%% Read that catalog
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include));
FT = T.FT(logical(T.include));
LT = T.LT(logical(T.include));
Kindex = find(T.include);


if ismember('FTa',T.Properties.VariableNames)
    FT = T.FTa(logical(T.include));
    LT = T.LTa(logical(T.include));
end

%% Concatenate the catalogued data [X,Y,n_bins]
% We need equal number of trials for each experiment

% Preallocate
traindata = cell(length(KWIKfiles),length(BinSizes));
n_bins = nan(size(BinSizes));

% Collect binned data, flattened into a column
for R = 1:length(KWIKfiles)
    efd = efds(R);
    Trials = FT(R):LT(R);
    
    if strcmp(T.VOI(Kindex(R)),'A')
        VOI = [2:11];
    end
    
    [Raster_S] = efd.ValveSpikes.RasterAlign(:,1,:); % no scootchin
    Raster_S = squeeze(Raster_S(:,1,:));
    
    %     align to first sniff after offset
    for V = 1:size(Raster_S,1)
        for TR = 1:size(Raster_S{V,1},1)-1
            scootchtime = efd.PREX(find(efd.PREX > efd.ValveTimes.FVSwitchTimesOff{V}(TR),1)) - efd.ValveTimes.PREXTimes{V}(TR);
            for U = 1:size(Raster_S,2)
                Raster_S{V,U}{TR} = Raster_S{V,U}{TR} - scootchtime;
            end
        end
    end
    
    for BS = 1:length(BinSizes)
        [Y,traindata{R,BS},n_bins(BS)] = BinRearranger_KDF(Raster_S([1,VOI],:,:),PST,BinSizes(BS),Trials);
        
    end
end

% Concatenate across experiments.
for BS = 1:length(BinSizes)
    X{BS} = cat(2,traindata{:,BS});
end

end