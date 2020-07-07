function [SpikeMat,PSTHt] = dpcaPrepper_cross(Catalog, efds, BinSize, PST)

%% Read that catalog
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include));
Kindex = find(T.include);

statelist = 'ak';

for state = 1:length(statelist)
    clear Spikes
    
    FT = T.(['FT',statelist(state)]);
    FT = FT(Kindex);
    LT = T.(['LT',statelist(state)]);
    LT = LT(Kindex);
    
    %% Concatenate the catalogued data [X,Y,n_bins]
    
    % Preallocate
    traindata = cell(length(KWIKfiles),length(BinSize));
    n_bins = nan(size(BinSize));
    % Collect binned data, flattened into a column
    for R = 1:length(KWIKfiles)
        efd = efds(R);
        Trials = FT(R):LT(R);
        
        if strcmp(T.VOI(Kindex(R)),'A')
            VOI = [4,7,8,12,15,16];
        elseif strcmp(T.VOI(Kindex(R)),'C')
            VOI = [11,7,8,6,12,10];
        end
        
        %%
        
        [~, PSTHtrials, PSTHt] = PSTHmaker_Beast(efd.ValveSpikes.RasterAlign([1,VOI],:,:), PST, BinSize, Trials);
        
        for N = 1:size(PSTHtrials,2)
            for S = 1:size(PSTHtrials,1)
                for t = 1:size(PSTHtrials,3)
                    Spikes{R}(N,S,:,t) = PSTHtrials{S,N,t};
                end
            end
        end
        
    end
    
    %% Concatentate across experiments
    SpikeMat{state} = cat(1,Spikes{:});
end

%% Do celltyping.
Type = 'Stable';
specialparams.FRlim = 1/100;
specialparams.UVlim = 50;
specialparams.DFRlim = 100;
[~, TypeStack] = CellTyper(Catalog, Type, specialparams);

for state = 1:2
    SpikeMat{state} = SpikeMat{state}(TypeStack{1},:,:,:,:);
end

% Concatenate across experiments.
% for BS = 1:length(BinSizes)
%     X{BS} = cat(2,traindata{:,BS});
% end

end