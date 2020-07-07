function [ypos,pttime,asym,hfw,bigwave,ypos_real,ptslp,ptrat_1,ptrat_2] = WaveStats_Beast(Wave)
pos = cell2mat(Wave.Position');
ypos = pos(:,2);
ypos_real = ypos;
ypos = ypos-median(ypos);

clear bigwave
x = Wave.AverageWaveform(:);
for k = 1:length(x)
    %     [~,b] = max(peak2peak(x{k}'));
    [~,b] = min(min(x{k}'));
    
    bigwv = x{k}(b,:);
    bigwave(k,:) = interp(bigwv,4);
end

for k = 1:size(bigwave,1)
    [tro,troloc] = min(bigwave(k,:));
    [pk2,pk2loc] = max(bigwave(k,troloc:end)+0);
%     pk2 = pk2-tro;
    [pk1,pk1loc] = max(bigwave(k,1:troloc)+0);
%     pk1 = pk1-tro;
    pttime(k) = (1/4)*(1/30)*pk2loc;
    asym(k) = (pk2-pk1)/(pk2+pk1);
    halfheight = 0.5*tro;
%         hfw(k) = (1/4)*(1/30)*(find(bigwave(k,troloc:end)>halfheight,1)+troloc-find(bigwave(k,1:troloc)<halfheight,1));


    ptrat_1(k) = abs(pk1)/abs(tro);
    ptrat_2(k) = abs(pk2)/abs(tro);
    ptslp(k) = (pk2-tro)/pttime(k);
    %     auc(k) = sum((bigwave(troloc:end)>0).*bigwave(troloc:end))/tro;
    
    
    % spikewidth
    Fs = 30000*4; % sampling freq
    wintype = 'hanning';    % spectrum
    nfft = 1024;            % spectrum
    nsamples = size(bigwave(k,:),2); % nb of samples per wavform
    win = hanning( nsamples );
    
    [S,F] = periodogram(bigwave(k,:),[],1024,Fs);
    [pk,loc] = max(S);
    hfw(k) = 1000/F(loc);
    
     % From spikes_stats (Eran):
        [ pow2 f ] = my_spectrum( bigwave(k,:)', nfft, Fs, win, 0, 0 );
        
        mp = mean( pow2( 2 : end, : ), 2 );
        [ ign, midx ] = max( mp );
        fmax = f( midx ); % peak of spectrum
        hfw(k) = 1/fmax*1000;
    
end
% 
% WVstats.ypos = ypos;
pttime = pttime';
asym = asym';
hfw = hfw';
ptslp = ptslp';
ptrat_1 = ptrat_1';
ptrat_2 = ptrat_2';
end
