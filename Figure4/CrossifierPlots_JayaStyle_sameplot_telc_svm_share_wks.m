clearvars
close all
clc

%% Parameters
% You will need to feed the function things.
% Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_pcx_awk_kx_telcontrol_stabletrials.txt';
Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_pcx_awk_kx_telc_stabletrials.txt';
% Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_bulb_awk_kx_telcontrol_stabletrials.txt';
% Catalog = 'Z:\expt_sets\catalogs\AK2\ExperimentCatalog_bulb_awk_kx_telc_stabletrials.txt';

BinSizes = 0.48;
PST = [0 0.485];

Permutations = 50; % should be 200

samplesizelist = [10:10:170]; % for pcx
% samplesizelist = [5:5:60]; % for bulb

specialparams.FRlim = 1/100;
specialparams.UVlim = 50;
specialparams.DFRlim = 100;

%% defining tasks
task = 'default';
Type = 'Stable';

%% Run Flexifier
[ACC, ACCx, TypeStack, task, p, X ,Y, SelectCells] = Crossifier_liblinear(Catalog,BinSizes,PST,Permutations,samplesizelist,task,Type,specialparams);

%% plotting
figure(1)
printpos([300 300 260 180])
clf
colors = [0 0 0; 0.051 0.447 0.7294];

for cty = 1:size(TypeStack,2)
    for tk = 1:length(task)
        
        %% within state
        subplotpos(2,1,1,1,.2)
        for state = 1:2
            for BS = 1:length(BinSizes)
                ssl = samplesizelist;
                ssl(isnan(nanmean(ACC{state,cty,BS,tk},2))) = nan;
                
                
                 %% bootstrap on means
                tosample = ACC{state,cty,BS,tk};
                clear bsmn
                for bts = 1:1000
                bssample = randi(size(tosample,2),size(tosample,2),1);
                bsmn(:,bts) = mean(tosample(:,bssample),2);
                end
                
                co = prctile(bsmn',[.5, 99.5]);
                
                [l,p] = boundedline(ssl,nanmean(ACC{state,cty,BS,tk},2),[nanmean(ACC{state,cty,BS,tk},2) - co(1,:)', co(2,:)'- nanmean(ACC{state,cty,BS,tk},2)],'cmap',colors(state,:),'alpha'); hold on;
                xlim([0 max(samplesizelist)+10])
                box off; axis square;
                set(gca,'Clipping','off')
                
            end
        end
        %%
        
        xx = ACC{1,cty,BS,tk}-ACC{2,cty,BS,tk};
        plot([0 max(samplesizelist)],[100/6 100/6],'k--') % zero line
        
        % axis cleanup
        ylim([0 100])
        xlim([0 max(samplesizelist)])
        box off; axis square;
        set(gca,'Clipping','off','YTick',[0 100],'YTickLabel',{'0', '100'},'XTick',0:100:400)
        yl = ylabel('accuracy (%)');
        yl.Position(1) = -40;
        xlabel('number of cells')
        set(gca,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1.2)
        %%
        
        % cross state
        subplotpos(2,1,2,1,.2)
        for state = 2
            for BS = 1:length(BinSizes)
                ssl = samplesizelist;
                ssl(isnan(nanmean(ACCx{state,cty,BS,tk},2))) = nan;
                                
                %% bootstrap on means
                tosample = ACCx{state,cty,BS,tk};
                clear bsmn
                for bts = 1:1000
                bssample = randi(size(tosample,2),size(tosample,2),1);
                bsmn(:,bts) = mean(tosample(:,bssample),2);
                end
                
                co = prctile(bsmn',[.5, 99.5]);
                [l,p] = boundedline(ssl,nanmean(ACCx{state,cty,BS,tk},2),[nanmean(ACCx{state,cty,BS,tk},2) - co(1,:)', co(2,:)'- nanmean(ACCx{state,cty,BS,tk},2)],'alpha','cmap',colors(state,:)); hold on;
                ylim([0 1]); set(gca,'YTick',[0 .5 1])
                xlim([0 max(samplesizelist)])
                box off; axis square;
                set(gca,'Clipping','off')
                
            end
        end
        
        xx = ACCx{1,cty,BS,tk}-ACCx{2,cty,BS,tk};
        co = prctile(xx',[2.5, 97.5]);
        
        plot([0 max(samplesizelist)],[100/6 100/6],'k--') % zero line
        
        % axis cleanup
        ylim([0 60])
        xlim([0 max(samplesizelist)])
        box off; axis square;
        set(gca,'Clipping','off','YTick',[0 60],'YTickLabel',[0 60],'XTick',0:100:400)
        xlabel('number of cells')
        set(gca,'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1.2)
    end
end