function FilesKK = FindFilesKK(KWIKfile)

FilesKK.KWIK = (KWIKfile);
[a,b] = fileparts(KWIKfile);
x = strfind(b,'.');
if ~isempty(x)
    FilesKK.KWX = [a,filesep,b(1:x(end)-1),'.kwx'];
    FilesKK.DAT = [a,filesep,b(1:x(end)-1),'.dat'];
    FilesKK.LFP = [a,filesep,b(1:x(end)-1),'.lfp'];
else
    FilesKK.KWX = [a,filesep,b,'.kwx'];
    FilesKK.DAT = [a,filesep,b,'.dat'];
    FilesKK.LFP = [a,filesep,b,'.lfp'];
end
x = strfind(b,'_');
FilesKK.AIP = [a,filesep,b(1:x(end)-1),'.ns3'];

[a,b] = fileparts(KWIKfile);
if length(strfind(KWIKfile,'.'))<2
    FilesKK.EFD = [a,filesep,b,'.efd'];
else
    FilesKK.EFD = [a,filesep,b(1:strfind(b,'.')),'efd'];
end