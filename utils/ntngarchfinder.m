function [AP, AN] = ntngarchfinder(KWIKfile)
[efd] = loadEFD(KWIKfile);
[SpikeTimes] = SpikeTimes_Beast(FindFilesKK(KWIKfile));

[ypos,pttime] = WaveStats_Beast(SpikeTimes.Wave);

trials = 1:20;
[~,PSTHtrials,PSTHt] = PSTHmaker_Beast(efd.LaserSpikes.RasterAlign,[-.5 .9],.001,trials);

ezptrials = cellcat3d(squeeze(PSTHtrials(1,:,:)));

for unit = 1:size(efd.LaserSpikes.SpikesBeforeLaser,2)
    Control = efd.LaserSpikes.SpikesBeforeLaser{unit};
    Stimulus = efd.LaserSpikes.SpikesDuringLaser{unit};
    [~, p(unit)] = RankSumROC(Control, Stimulus); % ranksum select laser responsive cells
    Rate(unit) = (mean(Control));
    
    %% median last spike
    temp = squeeze(ezptrials(unit,:,PSTHt>0));
    for tr = 1:size(temp,1)
        templat = find(temp(tr,:),1,'last');
        if isempty(templat)
            lat(tr) = 0;
        else
            lat(tr) = 0.0005*templat;
        end
    end
    medlat(unit) = median(lat);
end

disq = pttime'< 0.35 | Rate < .175;
AP = ~disq & medlat<0.01 & p<.0001;
AN = ~disq & (medlat>0.01 | p>.0001);
AD = disq;