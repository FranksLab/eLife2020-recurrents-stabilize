clearvars
close all
clc

%% Parameters
Catalog = 'ExperimentCatalog_ntng_kx_share.txt';

BinSize = 0.2;
PST = [-1.5 4];

% defining tasks
task = 'default';
Type = 'StableNTNG';

% for all units
specialparams.FRlim = 0;
specialparams.UVlim = 1000000;
specialparams.DFRlim = 10000000;

% WindowType = 'Expanding';
Permutations = 100;
samplesize = 100;

%% Run Flexifier
[ACCpf, TypeStack, task, win_t] = Flexifier_CrossTime_loadingSets_permute_loo_ntng(Catalog,BinSize,PST,Permutations,samplesize,'Sliding',task,Type,specialparams);
%% Postprocessing
figure(18)
printpos([200 200 900 300])
clf

colores = [1 0 0; 0 0 0];

win_t = linspace(PST(1),PST(2),size(ACCpf{1,1},1));

for ct = 1:2
    for p = 1:Permutations
        acc_fix(p,:,:) = ACCpf{ct,p};
        acc_diag(p,:) = diag(squeeze(acc_fix(p,:,:)));
    end
    
     cmap{6} = colorcet('R2');
    
    subplot(1,4,ct)
    c = 6;
    imagesc(win_t,win_t,squeeze(mean(acc_fix)))
    axis square
    caxis([1/6 1])
    set(gca,'XTick',[-1:1:4],'YTick',[-1:1:4])
    colormap(cmap{c});
    xlabel('trained time bin')
    ylabel('tested time bin')
    
    subplot(1,4,3)
    hold on
    if ct == 1
        area([-1.2 -.2], [1 1], 'edgecolor', 'none', 'facecolor', [.7 .9 1],'facealpha', .5)
    end
    
    m = bootstrp(1000,@mean,acc_diag);
    
    plotvar = squeeze(mean(acc_fix));
    % boundedline(win_t,mean(acc_diag),prctile(m,2.5)-mean(acc_diag),mean(acc_diag))
    
    co = prctile(m,[.5, 99.5]);
    eb = [mean(acc_diag) - co(1,:); co(2,:)- mean(acc_diag)]';
    [l,p] = boundedline(win_t,mean(acc_diag),eb,'cmap',colores(ct,:),'alpha'); hold on;
    hold on
    xlim([PST])
    ylim([0 1])
    axis square
    box off
    xlabel('time (s)')
    ylabel('accuracy')
    set(gca,'XTick',[-1:1:4])
    
    
    functiontimelength = 3;
     for p = 1:Permutations
         plotvar = squeeze(acc_fix(p,:,:));
         temp{ct} = tril(plotvar);
         temp{ct}(temp{ct} == 0) = nan;
         for k = 1:size(temp{ct},1)
             temp{ct}(:,k) = circshift(temp{ct}(:,k),-k+1);
         end
         
         temp{ct} = temp{ct}(win_t-PST(1)<=functiontimelength, win_t>=0 & PST(2)-win_t>=functiontimelength); %length of function [3], start points [0 to end of window - enough bins to have a full function]
         for k = 1:size(temp{ct},2)
             sum_decay(ct,k,p) = mean(temp{ct}(2:end,k));
         end
         
         decaycurve(ct,:,p) = temp{ct}(:,1);
     end
end

time = win_t(win_t>=0 & PST(2)-win_t>=functiontimelength);

subplot(1,4,4)
cla
for ct = 1:2
    delta_t = win_t(1:size(temp{ct}(:,k),1))-win_t(1);
    m = bootstrp(1000,@mean,squeeze(decaycurve(ct,:,:))');
    co = prctile(m,[.5, 99.5]);
    eb = [mean(squeeze(decaycurve(ct,:,:))') - co(1,:); co(2,:)- mean(squeeze(decaycurve(ct,:,:))')]';
    [l,p] = boundedline(delta_t,mean(squeeze(decaycurve(ct,:,:))'),eb,'cmap',colores(ct,:),'alpha'); hold on;
end
axis square
ylim([0 1])
box off
ylabel('accuracy')
xlabel('\Deltat')