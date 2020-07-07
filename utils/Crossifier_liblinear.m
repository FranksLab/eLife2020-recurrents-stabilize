function [ACC, ACCx, TypeStack, task, p, X, Y, SelectCells] = Crossifier_liblinear(Catalog,BinSizes,PST,Permutations,samplesizelist,task,Type,specialparams)

%% default task
if strcmp(task,'default')
    clear task
    %% Defining tasks
    % identity coding
    task{1}.taskstim{1} = [2:7];
    
end

%% Do celltyping.
[~, TypeStack] = CellTyper(Catalog, Type, specialparams);

%% Load in all EFD data to save on processing time.
efds = EFDloader(Catalog);

%% Concatenate the catalogued data [X,Y,n_bins] = PseudoPopulator(Catalog, BinSizes, PST)
statelist = ['a','k'];
[X,Y,n_bins] = CrossPopulator(Catalog, efds, BinSizes, PST, statelist); % if abs(BinSize-diff(PST))<BinSize there will be 1 bin.

if isfield(specialparams,'ZBin')
    for state = 1:size(X,1)
        for BS = 1:size(X,2)
            temp =  X{state,BS}(Y{state}==1,:);
            threshtemp = mean(temp)+std(temp);
            X{state,BS} = bsxfun(@gt, X{state,BS},threshtemp);
        end
    end
end

if isfield(specialparams,'Z')
    for BS = 1:size(X,2)
        for state = 1:size(X,1)
            temp =  X{state,BS}(Y{state}==1,:);
            stdtemp{state} = std(temp);
            meantemp{state} = mean(temp);
            
            %%
            X{state,BS} = bsxfun(@minus, X{state,BS}, meantemp{state});
            X{state,BS} = bsxfun(@rdivide, X{state,BS}, stdtemp{state});
            
            X{state,BS}(:,stdtemp{state}==0) = 0;
        end
    end
end

if isfield(specialparams,'mnsubtract')
    for BS = 1:size(X,2)
        silents = [];
        for state = 1:size(X,1)
            temp =  X{state,BS}(Y{state}==1,:);
            stdtemp{state} = std(temp);
            meantemp{state} = mean(temp);
            
            %%
            X{state,BS} = bsxfun(@minus, X{state,BS}, meantemp{state});            
        end
    end
end

%% Multiple permutations
for pm = 1:Permutations
    
    %% Loop through sample sizes and do subselection as a function, then do
    % classifying as a function for each sample size
    for ssz = 1:length(samplesizelist)
        samplesize = samplesizelist(ssz);
        % Subselect cells/bins (TypeStack, BinSizes, n_bins)
        [Xs, SelectCells{ssz,pm}] = SubsetPermuter_Cross(X, Y, samplesize, TypeStack, BinSizes, n_bins);
        
        for cty = 1:size(TypeStack,2)
            for BS = 1:length(BinSizes)
                for tk = 1:length(task)
                    for state = 1:2
                        if sum(cellfun(@isempty,Xs(cty,BS,state)))>0
                            ACCpf{pm}(state,cty,BS,tk,ssz) = nan;
                            Npf{pm}(state,cty,BS,tk,ssz) = nan;
                            
                        else
                            [ACCpf{pm}(state,cty,BS,tk,ssz), model{pm}{state}] = GenTaskClassifier_liblinear(Y{state}, Xs{cty,BS,state}, task{tk});
                            Npf{pm}(state,cty,BS,tk,ssz) = nan;
                        end
                    end
%                     ACCx = nan;
                    for state = 1:2
                        if sum(cellfun(@isempty,Xs(cty,BS,state)))>0
                            ACCpfX{pm}(state,cty,BS,tk,ssz) = nan;
                            NpfX{pm}(state,cty,BS,tk,ssz) = nan;
                        else
                            A = model{pm}{(1:2)~=state};
                            testidx = ismember(Y{state},task{tk}.taskstim{1});
                            
                            [~, accuracy, ~] = predict(Y{state}(testidx), sparse(Xs{cty,BS,state}(testidx,:)), A);

                            ACCpfX{pm}(state,cty,BS,tk,ssz) = accuracy(1);
                            NpfX{pm}(state,cty,BS,tk,ssz) = nan;
                        end
                    end
                end
            end
        end
    end
end

%% Rearrange the awkward data from the parfor loop
ACCpf = cat(6,ACCpf{:});
Npf = cat(6,Npf{:});
for state = 1:2
    for cty = 1:size(TypeStack,2)
        for BS = 1:length(BinSizes)
            for tk = 1:length(task)
                ACC{state,cty,BS,tk} = squeeze(ACCpf(state,cty,BS,tk,:,:));
                N{state,cty,BS,tk} = squeeze(Npf(state,cty,BS,tk,:,:));
            end
        end
    end
end

ACCpfX = cat(6,ACCpfX{:});
NpfX = cat(6,NpfX{:});
for state = 1:2
    for cty = 1:size(TypeStack,2)
        for BS = 1:length(BinSizes)
            for tk = 1:length(task)
                ACCx{state,cty,BS,tk} = squeeze(ACCpfX(state,cty,BS,tk,:,:));
                Nx{state,cty,BS,tk} = squeeze(NpfX(state,cty,BS,tk,:,:));
            end
        end
    end
end

p = nan;
end