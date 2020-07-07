function [AP, AN] = ntngchr2finder(KWIKfile, pfourfive)
[efd] = loadEFD(KWIKfile);
[SpikeTimes] = SpikeTimes_Beast(FindFilesKK(KWIKfile));

[~,pttime] = WaveStats_Beast(SpikeTimes.Wave);

LBlist = [1 find(abs(diff(efd.LaserTimes.LaserOn{1},2))>.01)+2 length(efd.LaserSpikes.RasterAlign{1})+1];
LBlist(diff(LBlist)==1) = [];
clear PSTH PSTHtrials ezptrials

% trials = LBlist(pfourfive(1)):(LBlist(pfourfive+1+1)-1);
[PSTH,PSTHtrials,PSTHt] = PSTHmaker_Beast(efd.LaserSpikes.RasterAlign,[0 .01],.0005,LBlist(pfourfive(1)):(LBlist(pfourfive(1)+1+1)-1));

ezptrials = cellcat3d(squeeze(PSTHtrials(1,:,:)));
test = cellcat3d(PSTH);


[pk,lc] = max(squeeze(test)');
latency = PSTHt(lc);
latency(lc==1) = nan;

clear prePSTHtrials
[~, prelaseridx] = max(diff(efd.LaserTimes.LaserOn{1}));
%
clear P
for unit = 1:length(SpikeTimes.tsec)
    [~,P{unit}] = PSTHmaker_Beast({{SpikeTimes.tsec{unit}}},[efd.LaserTimes.LaserOn{1}(prelaseridx+1)-60 efd.LaserTimes.LaserOn{1}(prelaseridx+1)-0],.0005,1);
end
PP = cell2mat(cat(1,P{:})');

% SALT
clear pSALT ISALT
for unit = 1:length(SpikeTimes.tsec)
    
    spt_baseline = reshape(PP(:,unit),[],100)';
    spt_test = squeeze(ezptrials(unit,:,:));
    
    [pSALT(unit) ~] = salt(spt_baseline,spt_test,0.0005);
end

disq = pttime'< 0.35;

AP = ~disq & latency <= 0.003 & pSALT<=.001;
AN = ~disq & (latency > 0.003 | pSALT>.001);
AD = disq;