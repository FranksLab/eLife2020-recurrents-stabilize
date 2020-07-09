function [CM, ACC] = SubsetCrossTimeClassifier_liblinear(trainlabel, traindata, task)

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
    tempdata = tempdata(ismember(trainlabel,task.taskstim{j}),:,:);
    
    %     tempdata = shake(tempdata,1);
    %     templabel = shake(templabel,1);
    obsindex = 1:length(templabel);
    
    for o = obsindex
        for bin = 1:size(tempdata,2)
            
            trl = templabel(obsindex ~= o);
            trd = tempdata(obsindex~=o,bin,:);
            clslist = unique(trl,'stable');
            
            
            %             ACC = train(trl,sparse(squeeze(trd)),'-c 1 -q -s 4 -v 10');
            model = train(trl,sparse(squeeze(trd)),'-c 1 -q -s 4');
            
            
            for binX = 1:size(tempdata,2)
                ted = tempdata(o,binX,:);
                tel = templabel(o);
                
                %                 A = model{pm}{(1:2)~=state};
                %                 testidx = ismember(Y{state},task{tk}.taskstim{1});
                %
                %                 [pred(o,bin,binX), ~, ~] = predict(tel, sparse(squeeze(ted)), model);
                [pred(o,bin,binX), ~, ~] = predict(tel, sparse(squeeze(ted)'), model);
                
                %
                %                 % testing
                %                 distances = pdist([squeeze(tempdata(o,binX,:))';clsmean],DISTANCE);
                %                 distances = distances(1:length(clslist));
                %                 if isfield(task,'notes')
                %                     if strcmp(task.notes,'NRA') && ismember(find(clslist==templabel(o)),task.fixup{j})
                %                         distances(clslist == templabel(o)) = nan;
                %                     end
                %                 end
                %                 [~, pred(o,bin,binX)] = min(distances);
            end
        end
    end
    
    %     pred = clslist(pred);
    
    % calculate accuracy
    for bin = 1:size(tempdata,2)
        for binX = 1:size(tempdata,2)
            CM{j}(bin,binX,:,:) = confusionmat(templabel,squeeze(pred(:,bin,binX)));
            CM{j}(bin,binX,:,:) = bsxfun(@rdivide,squeeze(CM{j}(bin,binX,:,:)),sum(squeeze(CM{j}(bin,binX,:,:)),2));
            ACC(bin,binX) = mean(diag(squeeze(CM{j}(bin,binX,:,:))));
        end
    end
    %
    %     % calculate accuracy. This will become more complicated
    %         CM{j} = confusionmat(templabel,pred);
    %         CM{j} = bsxfun(@rdivide,CM{j},sum(CM{j},2));
    %         ACC(j) = mean(diag(CM{j}));
    %         N(j) = length(templabel);
    %
end

% ACC = mean(ACC);

end

