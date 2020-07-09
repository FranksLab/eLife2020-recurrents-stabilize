clearvars
close all
clc

specialparams.FRlim = 1/100;
specialparams.UVlim = 50;
specialparams.DFRlim = 100;

%% Pick out files with 'kwik' in its name and put each in one cell
Catalog = 'ExperimentCatalog_bulb_awk_kx_share.txt';
T = readtable(Catalog, 'Delimiter', ' ');
Ofiles = T.kwikfile(logical(T.include));
[TypeIdxO, ~] = CellTyper (Catalog, 'Stable', specialparams);


Catalog = 'ExperimentCatalog_pcx_awk_kx_share.txt';
T = readtable(Catalog, 'Delimiter', ' ');
Pfiles = T.kwikfile(logical(T.include));
[TypeIdxP, ~] = CellTyper (Catalog, 'Stable', specialparams);

%% Choose an example
exa = 4;

%% Get respiration info
FilesKK = FindFilesKK(Pfiles{exa});
[a,b] = fileparts(FilesKK.AIP);
RESPfile = [a,filesep,b,'.resp'];
[InhTimes,PREX,POSTX,RRR,BbyB] = resp_loader(RESPfile);

[Fs,t,VLOs,FVO,resp,LASER,PID] = NS3Unpacker(FilesKK.AIP);
% RRR = PID;

%% Choose LFP data channels
% Bulb
[SpikeTimes_O] = SpikeTimesPro(FindFilesKK(Ofiles{exa}));
[channelposlist, nchannels] = ChannelPositionExtractor(Ofiles{exa});

% Find channel closest to mean of all waveforms positions.
[WFxy] = mean(cell2mat(SpikeTimes_O.Wave.Position(2:end)'));
D = pdist2(WFxy,channelposlist);
[~, SomaChan_O] = min(D);

% Pcx
[SpikeTimes_P] = SpikeTimesPro(FindFilesKK(Pfiles{exa}));
[channelposlist, nchannels] = ChannelPositionExtractor(Pfiles{exa});

% Find channel closest to mean of all waveforms positions.
[WFxy] = mean(cell2mat(SpikeTimes_P.Wave.Position(2:end)'));
D = pdist2(WFxy,channelposlist);
[~, SomaChan_P] = min(D);

%% Find times to censor near odor presentations
efd = loadEFD(Pfiles{exa});
Raster = efd.ValveSpikes.RasterAlign;
FVall = cat(2,efd.ValveTimes.FVSwitchTimesOn{:});
FVall = sort(FVall);
SafeWindowStarts = FVall+4;
OdorWindowStarts = FVall;

%% Build bulb and pcx population rasters.
[SpikeTimes_O] = SpikeTimesPro(FindFilesKK(Ofiles{exa}));
pos = cell2mat(SpikeTimes_O.Wave.Position(2:end)');
[sortpos, posidx] = sort(pos(:,2));
clear b
for k = 1:size(SpikeTimes_O.tsec,1)
    b{k} = SpikeTimes_O.tsec{k}';
end
PopRaster_O = b(2:end);

[SpikeTimes_P] = SpikeTimesPro(FindFilesKK(Pfiles{exa}));
pos = cell2mat(SpikeTimes_P.Wave.Position(2:end)');
[sortpos, posidx] = sort(pos(:,2));
clear b
for k = 1:size(SpikeTimes_P.tsec,1)
    b{k} = SpikeTimes_P.tsec{k}';
end
PopRaster_P = b(2:end);

%%
figure(1)
set(0,'defaultlinelinewidth',1.0)
set(0,'defaultaxeslinewidth',1.0)
printpos([100 100 1000 250]);
clf

mm = 28;
exmin = [(OdorWindowStarts(197)), (OdorWindowStarts(337)), (OdorWindowStarts(834))];

breathlim = [-4000 4000];

for snip = 1:2

% Examples
exwin = [exmin(snip)-2 exmin(snip)+2.5];
subplotpos(5,3,1,snip,.2)
plot(exwin(1):1/2000:exwin(2),RRR(ceil(exwin(1)*2000:exwin(2)*2000)),'k')
xlim(exwin)
set(gca,'XTick',[],'YTick',[])
hold on;

subplotpos(5,3,2,snip,.2)
% read in bulb LFP data
[datmat] = lfp_reader(Ofiles{exa},0,exwin,nchannels);
datmat = datmat(SomaChan_O,:)/4000; % now it's in millivolts
plot(datmat,'k')
ylim([-1 1])
xlim([0 length(datmat)])
set(gca,'XTick',[],'YTick',[])
box on

subplotpos(5,3,3,snip,.2)
% read in pcx LFP data
[datmat] = lfp_reader(Pfiles{exa},0,exwin,nchannels);
datmat = datmat(SomaChan_P,:)/4000; % now it's in millivolts
plot(datmat,'k')
ylim([-1 1])
xlim([0 length(datmat)])
set(gca,'XTick',[],'YTick',[])
box on

fibo = find(TypeIdxO{exa}(2:end));
fibp = find(TypeIdxP{exa}(2:end));


subplotpos(5,3,4,snip,.2)
plotSpikeRaster(PopRaster_O(fibo), 'PlotType','vertline','XLimForCell',exwin,'VertSpikeHeight',.7);
axis xy
ylim([0 length(PopRaster_O(fibo))+1])
set(gca,'XTick',[],'YTick',[])
box on
hold on
plot([exwin(1)+2, exwin(1)+2],[1 30],'r')

% 
subplotpos(5,3,5,snip,.2)
plotSpikeRaster(PopRaster_P(fibp), 'PlotType','vertline','XLimForCell',exwin,'VertSpikeHeight',.7);
axis xy
ylim([0 length(PopRaster_P(fibp))+1])
set(gca,'XTick',[],'YTick',[])
box on
hold on
plot([exwin(1)+2, exwin(1)+2],[1 70],'r')

end

