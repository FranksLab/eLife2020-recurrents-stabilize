function [Scores,efd] = SCOmaker_NoBlank(KWIKfile,TrialSets)
% function [Scores,efd] = SCOmaker(KWIKfile,TrialSets)
% This SCOmaker function will give you information for all odor-cell pairs
% about how the odor reponse differs from the control response. The output
% will be a structure called Scores. These can be further summarized to
% characterize unit responses per experiment. It will also pass along the
% full experiment data (efd) gathered from EFDmaker.
%
% Scores output are stored like this:
% Scores.BlankRate and Scores.BlankSD - (Unit,BreathCycle,TrialSet)
% All other Scores - (Valve,Unit,BreathCycle,TrialSEt)
% BreathCycle refers to which respiration cycle after odor delivery is
% analyzed. Usually you want to analyze BreathCycle 1 but you could
% potentially care about the evolution of a response past the first cycle.

%% Here we are gathering information.
% Basically spike counts relative to valve openings and when experimental events occur.
[efd] = loadEFD(KWIKfile);

% Also if the user doesn't specify any trials to analyze in particular we
% will find the minimum number of trials for a single valve and use those.
mintrials = min(cellfun(@length,efd.ValveTimes.FVSwitchTimesOn));
if nargin < 2
    TrialSets{1} = 1:mintrials;
end

%% Here we will take any measure of the response and give indications of
% how different the odorant response is from valve 1.
% These will include: auROC, p-value for ranksum test, z-score, rate
% change, mean and SD of valve 1 response.

%% Cycles
TESTVAR = efd.ValveSpikes.MultiCycleSpikeRate;
[SpikeTimes] = SpikeTimesPro(FindFilesKK(KWIKfile));
[~,BLANKVAR,~] = VSMultiCycleCount(efd.ValveTimes,SpikeTimes,efd.PREX,{-2});

for tset = 1:length(TrialSets)
    if ~isempty(TrialSets{tset})
        for Unit = 1:size(TESTVAR,2)
            for Cycle = 1:size(TESTVAR,3)
                
                for Valve = 1:size(TESTVAR,1)
                    if ~isempty(efd.ValveTimes.PREXTimes{Valve})
                        
                        % Blank Rate and SD
                        Scores.BlankRate(Valve,Unit,Cycle,tset) = nanmean(BLANKVAR{Valve,Unit,1}(TrialSets{tset}));
                        Scores.BlankSD(Valve,Unit,Cycle,tset) = nanstd(BLANKVAR{Valve,Unit,1}(TrialSets{tset}));
                        
                        % First Cycle
                        % auROC and p-value for ranksum test
                        [Scores.auROC(Valve,Unit,Cycle,tset), Scores.AURp(Valve,Unit,Cycle,tset)] = RankSumROC(BLANKVAR{Valve,Unit,1}(TrialSets{tset}),TESTVAR{Valve,Unit,Cycle}(TrialSets{tset}));
                        
                        % Z-Scores based on valve 1 responses vs everything else.
                        Scores.ZScore(Valve,Unit,Cycle,tset) = (nanmean(TESTVAR{Valve,Unit,Cycle}(TrialSets{tset}))-nanmean(BLANKVAR{Valve,Unit,1}(TrialSets{tset})))./nanstd(BLANKVAR{Valve,Unit,1}(TrialSets{tset}));
                        Scores.ZScore(isinf(Scores.ZScore)) = NaN;
                        
                        % Rate change based on valve 1 responses vs everything else. The
                        % denominator comes into play here.
                        Scores.RateChange(Valve,Unit,Cycle,tset) = (nanmean(TESTVAR{Valve,Unit,Cycle}(TrialSets{tset}))-nanmean(BLANKVAR{Valve,Unit,1}(TrialSets{tset})));
                        
                        % Raw Rate
                        Scores.RawRate(Valve,Unit,Cycle,tset) = nanmean(TESTVAR{Valve,Unit,Cycle}(TrialSets{tset}));
                        
                        % Fano Factor - variability compared to Poisson
                        Scores.Fano(Valve,Unit,Cycle,tset) = nanvar(TESTVAR{Valve,Unit,Cycle}(TrialSets{tset}))./nanmean(TESTVAR{Valve,Unit,Cycle}(TrialSets{tset}));
                    end
                end
            end
        end
    end
end
end
