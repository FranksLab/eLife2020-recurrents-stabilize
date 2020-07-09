clearvars
close all
clc

% Catalog = 'ExperimentCatalog_bulb_awk_kx_share.txt';
Catalog = 'ExperimentCatalog_pcx_awk_kx_share.txt';

BinSize = [0.5];
PST = [0 .5];

%% Load in all EFD data to save on processing time.
efds = EFDloader(Catalog);

%%
[SpikeMat,time] = dpcaPrepper_cross(Catalog, efds, BinSize, PST);

%% Rearrange the SpikeMat a little more
Ss = [0 1 2 3 4 5 6];

for state = 1:2
    for T = 1:length(time)
        for S = 1:6
            V = find(Ss == S);
            if ~isempty(V)
                firingRates(:,S,state,T,:) = SpikeMat{state}(:,V,T,:);
            else
                firingRates(:,S,state,T,:) = nan(size(SpikeMat{state}(:,1,T,:)));
            end
        end
    end
end

% use all stimuli
firingRates = firingRates(:,:,:,:,:);
trialNum = sum(squeeze(~isnan(firingRates(:,:,:,1,:))),4);

%% moving into copy-paste from dpca_demo.m
% computing PSTHs
firingRatesAverage = nanmean(firingRates, 5);

%% Define parameter grouping

% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; herewe ignore the 1st dimension
% (neurons), i.e. we have the following parameters:
%    1 - stimulus
%    2 - state
combinedParams = {{1}, {2}, {[1 2]}};

margNames = {'Stimulus', 'State', 'S/State Interaction'};

margColours = [23 100 171; 187 20 25; 114 97 171]/256;

timeEvents = 0;

% check consistency between trialNum and firingRates
for n = 1:size(firingRates,1)
    for s = 1:size(firingRates,2)
        for d = 1:size(firingRates,3)
            assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
        end
    end
end

%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations
% in a .mat file with a given name. Once computed, you can simply load
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

ifSimultaneousRecording = false;

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, 8, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

%%
Xfull = firingRates;

X = Xfull(:,:)';
Xcen = bsxfun(@minus, X, mean(X));
dataDim = size(Xfull);
Z = Xcen * W;

stimC1 = find(whichMarg == 1,3);
stateC1 = find(whichMarg == 2,1);
componentsToPlot = [stimC1,stateC1];

Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]);

figure(99)
printpos([300 300 700 250])
clf
for S = 1:6
    for C = 1:2
        
        %% individual state plots
        withinax = subplot(1,3,C);
                edgeCT = [0 0 0; .4 .8 1]; % for white background
        
                CT = [parula(6)]; % for white background
        
        
        color = CT(S,:);
        edgecolor = edgeCT(C,:);
        hold on
        covarianceEllipse3D(mean(squeeze(Zfull(1:3,S,C,1,:))')',cov(squeeze(Zfull(1:3,S,C,1,:))'),S*.95);

        view(3)
        grid on
        axis square
        
        % ob-pcx
        xlim([-23 18])
        ylim([-15 15])
        zlim([-15 15])
        
        colormap(withinax,CT)
        xlabel('stimpc1')
        ylabel('stimpc2')
        zlabel('stimpc3')
        
        %% overlay plot
        ax = subplot(1,3,3);
        
        hold on
        covarianceEllipse3D(mean(squeeze(Zfull(1:3,S,C,1,:))')',cov(squeeze(Zfull(1:3,S,C,1,:))'),C);
        colormap(ax,edgeCT)
        
        hold on
        view(3)
        axis square
        
        % ob-pcx
        xlim([-23 18])
        ylim([-15 15])
        zlim([-15 15])

        xlabel('stimpc1')
        ylabel('stimpc2')
        zlabel('stimpc3')
        alpha(.4)
        grid on
        
    end
end
 
%% Matusita overlap
for S = 1:6
    C = 1;
    M0 = mean(squeeze(Zfull(1:3,S,C,1,:))')';
    S0 = cov(squeeze(Zfull(1:3,S,C,1,:))');
    
    C = 2;
    M1 = mean(squeeze(Zfull(1:3,S,C,1,:))')';
    S1 = cov(squeeze(Zfull(1:3,S,C,1,:))');
    
    % Matusita's overlap
    p = 3; % n dimensions
    Q = (2^(p/2) * det(S0)^.25 * det(S1)^.25) / det(S0+S1)^.5;
    R = exp(-.25*(M0-M1)'*inv(S0+S1)*(M0-M1));
    z(S) = Q*R;
end

%% Matusita's overlap plot - box and whisker for NN
OBnsF = [0.2962    0.5227    0.4090    0.4812    0.5611    0.5450];
PCXnsF = [0.5025    0.6820    0.8009    0.6086    0.8982    0.8149];

figure(44)
printpos([300 300 150 200])
clf
colores = [1 .2 .2; 0.4 0.4 0.4];

figure_boxplot([OBnsF;PCXnsF]',{[];[]},[],[],'horizontal',colores);
hold on
plotSpread([OBnsF;PCXnsF]','distributionColors',colores)
xlim([0 3])
set(gca,'Clipping','off','YTick',[0 1],'XTick',[])
ylim([0 1])
box off

