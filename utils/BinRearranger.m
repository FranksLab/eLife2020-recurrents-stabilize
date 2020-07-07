function [trainlabel,traindata,n_bins] = BinRearranger(Raster,PST,BinSize,Trials)

[~, PSTHtrials, PSTHt] = PSTHmaker_Beast(Raster, PST, BinSize, Trials);

%%
    A = cell2mat(PSTHtrials);

if ndims(A) == 3
    B = permute(A,[3,1,2]);
   
    traindata = reshape(B,size(B,1)*size(B,2),[]);
    trainlabel = repmat(1:size(A,1),size(PSTHtrials,3),1);
    trainlabel = trainlabel(:);
    
    n_bins = length(PSTHt);
else
    B = permute(A,[4,1,2,3]);
    
    traindata = reshape(B,size(B,1)*size(B,2),[]);
    trainlabel = repmat(1:size(A,1),size(PSTHtrials,4),1);
    trainlabel = trainlabel(:);
    
    n_bins = length(PSTHt);
end
end

