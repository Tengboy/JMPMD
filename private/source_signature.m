function [signature,files]=source_signature(root)
% Hash MATLAB source only. Generated results do not alter a run's signature.
paths=[dir(fullfile(root,'*.m'));dir(fullfile(root,'src','*.m')); ...
    dir(fullfile(root,'src','functions','*.m'));dir(fullfile(root,'private','*.m'))];
files=struct('file',{},'sha256',{});combined='';
for k=1:numel(paths)
    p=fullfile(paths(k).folder,paths(k).name);
    fid=fopen(p,'rb');assert(fid>=0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
    hash=hashbytes(bytes);rel=strrep(p(numel(root)+2:end),'\','/');
    files(k)=struct('file',rel,'sha256',hash);combined=[combined rel ':' hash newline]; %#ok<AGROW>
end
signature=hashbytes(unicode2native(combined,'UTF-8'));
end
function h=hashbytes(bytes)
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
