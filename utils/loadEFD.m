function [efd] = loadEFD(KWIKfile)

[a,b] = fileparts(KWIKfile);
if length(strfind(KWIKfile,'.'))<2
    EFDfile = [a,filesep,b,'.efd'];
else
    EFDfile = [a,filesep,b(1:strfind(b,'.')),'efd'];
end

load(EFDfile,'-mat')
