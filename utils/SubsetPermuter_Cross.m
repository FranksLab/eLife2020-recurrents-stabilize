function [Xs, SelectCells] = SubsetPermuter_Cross(X, Y, samplesize, TypeStack, BinSizes, n_bins)

% Preallocate
SelectBins = cell(size(TypeStack,2),length(BinSizes));
Xs = cell(size(TypeStack,2),length(BinSizes),size(X,1));
SelectCells = nan;

for cty = 1:size(TypeStack,2)
    TypeTemp = find(TypeStack{cty});

    if length(TypeTemp)>=samplesize
        SelectCells = TypeTemp(randperm(length(TypeTemp),samplesize));

        for BS = 1:length(BinSizes)
            SelectBins{cty,BS} = bsxfun(@plus,(SelectCells'-1).*n_bins(BS),(1:n_bins(BS))');
            SelectBins{cty,BS} = SelectBins{cty,BS}(:);
            for state = 1:size(X,1)
                Xs{cty,BS,state} = X{state,BS}(:,SelectBins{cty,BS});


            end
        end
    end
end
end