function [ACCpf, TypeStack, task, win_t] = Flexifier_CrossTime_loadingSets_permute_loo_ntng(Catalog,BinSize,PST,Permutations,samplesize,~,task,Type,specialparams)

%%

%% default task
if strcmp(task,'default')
    clear task
    %% Defining tasks
    % identity coding
    task{1}.taskstim{1} = [2 3 4 5 7 11];
end

%% Do celltyping.
[~, TypeStack] = CellTyper(Catalog, Type, specialparams);

%% Load in all EFD data to save on processing time.
efds = EFDloader(Catalog);

%% Loop through bin sizes and time limits
StepSize = BinSize/2;
windowend = PST(1):StepSize:PST(2);%.5-BinSize/2; % sliding window
win_t = windowend+BinSize/2;
win_t = win_t(1:end-1);

%% Multiple permutations

%% Concatenate the catalogued data [X,Y,n_bins] = PseudoPopulator(Catalog, BinSizes, PST)
[X,Y,n_bins] = PseudoPopulator_loadingSets_ntng(Catalog, efds, BinSize, PST);
XX = reshape(X{1},length(Y),n_bins,[]);
%%
parfor pm = 1:Permutations
    for ct = 1:2
        
        % Subselect cells/bins (TypeStack, BinSizes, n_bins)
        TypeTemp = find(TypeStack{ct});
        SelectCells = TypeTemp(randperm(length(TypeTemp),samplesize));
        Xs = XX(:,:,SelectCells);
        
        %Classify
        [CM, ACCpf{ct,pm}] = SubsetCrossTimeClassifier_liblinear(Y, Xs, task{1});
    end
end

end