function [Scores,efd] = SCOmaker_Beast_NoBlank(KWIKfile,TrialSets)

%% Here we are gathering information.
% Basically spike counts relative to valve openings and when experimental events occur.
[efd] = loadEFD(KWIKfile);

% Also if the user doesn't specify any trials to analyze in particular we
% will find the minimum number of trials for a single valve and use those.
mintrials = min(cellfun(@length,efd.ValveTimes.FVSwitchTimesOn));
if nargin < 2
    TrialSets{1} = 1:mintrials;
end

%% Cycles
TESTVAR = efd.ValveSpikes.MultiCycleSpikeRate(:,:,:,:);
%     TESTVAR_select = (TESTVAR(5:12,2:size(TESTVAR,2),:));

TESTVAR = efd.ValveSpikes.MultiCycleSpikeRate;
[SpikeTimes] = SpikeTimesPro(FindFilesKK(KWIKfile));
[~,BLANKVAR,~] = VSMultiCycleCount_Beast(efd.ValveTimes,SpikeTimes,efd.PREX,{-2});


for tset = 1:length(TrialSets)
    if ~isempty(TrialSets{tset})
        for Unit = 1:size(TESTVAR,3)
            for Conc = 1:size(TESTVAR,2)
                for Cycle = 1:size(TESTVAR,4)
                    % Blank Rate and SD
                    if ~isempty(efd.ValveTimes.PREXTimes{1,Conc})
                        Scores.BlankRate(Conc,Unit,Cycle,tset) = nanmean(TESTVAR{1,Conc,Unit,Cycle}(TrialSets{tset}));
                        Scores.BlankSD(Conc,Unit,Cycle,tset) = nanstd(TESTVAR{1,Conc,Unit,Cycle}(TrialSets{tset}));
                    end
                    for Valve = 1:size(TESTVAR,1)
                        if ~isempty(efd.ValveTimes.PREXTimes{Valve,Conc})
                            % First Cycle
                            [Scores.auROC(Valve,Conc,Unit,Cycle,tset), Scores.AURp(Valve,Conc,Unit,Cycle,tset)] = RankSumROC(BLANKVAR{Valve,Conc,Unit,1}(TrialSets{tset}),TESTVAR{Valve,Conc,Unit,Cycle}(TrialSets{tset}));

                            % Z-Scores based on valve 1 responses vs everything else.
                            Scores.ZScore(Valve,Conc,Unit,Cycle,tset) = (nanmean(TESTVAR{Valve,Conc,Unit,Cycle}(TrialSets{tset}))-nanmean(TESTVAR{1,end,Unit,Cycle}(TrialSets{tset})))./nanstd(TESTVAR{1,end,Unit,Cycle}(TrialSets{tset}));
                            Scores.ZScore(isinf(Scores.ZScore)) = NaN;
                            
                            % Rate change based on valve 1 responses vs everything else. The
                            % denominator comes into play here.
                            Scores.RateChange(Valve,Conc,Unit,Cycle,tset) = (nanmean(TESTVAR{Valve,Conc,Unit,Cycle}(TrialSets{tset}))-nanmean(TESTVAR{1,end,Unit,Cycle}(TrialSets{tset})));
                            
                            % Raw Rate
                            Scores.RawRate(Valve,Conc,Unit,Cycle,tset) = nanmean(TESTVAR{Valve,Conc,Unit,Cycle}(TrialSets{tset}));
                            
                            % Fano Factor - variability compared to Poisson
                            Scores.Fano(Valve,Conc,Unit,Cycle,tset) = nanvar(TESTVAR{Valve,Conc,Unit,Cycle}(TrialSets{tset}))./nanmean(TESTVAR{Valve,Conc,Unit,Cycle}(TrialSets{tset}));
                        end
                    end
                end
            end
        end
    end
end
end
