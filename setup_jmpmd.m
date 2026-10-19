function root=setup_jmpmd
%SETUP_JMPMD Add this release's numerical code; no external project needed.
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'src'),fullfile(root,'src','functions'),'-begin');
end
