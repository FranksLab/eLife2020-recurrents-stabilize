function [ACC, model] = GenTaskClassifier_liblinear(trainlabel, traindata, task)

if nargin<3
    DISTANCE = 'euclidean';
elseif ~isfield(task,'DISTANCE')
    DISTANCE = 'euclidean';
end

if isfield(task,'taskstim')
    num_jobs = length(task.taskstim);
else
    num_jobs = 1;
end

%%
for j = 1:num_jobs
    templabel = trainlabel;
    tempdata = traindata;
    
    templabel = templabel(ismember(trainlabel,task.taskstim{j}),:);
    tempdata = tempdata(ismember(trainlabel,task.taskstim{j}),:);
    
    ACC = train(templabel,sparse(tempdata),'-c 1 -q -s 4 -v 10');
    model = train(templabel,sparse(tempdata),'-c 1 -q -s 4');
end

end

