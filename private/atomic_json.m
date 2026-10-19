function atomic_json(path,value)
tmp=[path '.pending'];fid=fopen(tmp,'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value));clear cleanup
[ok,msg]=movefile(tmp,path,'f');assert(ok,msg);
end
