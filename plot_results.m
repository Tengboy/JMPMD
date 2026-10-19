function audit=plot_results(dataDir,outDir,varargin)
%PLOT_RESULTS Redraw main figures from published CSV or a completed new run.
% plot_results                         % bundled published results
% plot_results('runs/full','figures')   % your own newly simulated CSV
% Optional name/value Font: Chinese-capable font available on your machine.
root=setup_jmpmd;
if nargin<1||isempty(dataDir),dataDir=fullfile(root,'results','published');end
if nargin<2||isempty(outDir),outDir=fullfile(root,'replotted');end
p=inputParser;p.addParameter('Font','');p.parse(varargin{:});fontName=char(p.Results.Font);
if isempty(fontName)
    fonts=listfonts;choices={'SimSun','宋体','Noto Serif CJK SC','Noto Serif SC','Source Han Serif SC', ...
        'Microsoft YaHei','微软雅黑','Microsoft YaHei UI','PingFang SC','WenQuanYi Zen Hei','Noto Sans SC'};
    found=find(ismember(choices,fonts),1);
    assert(~isempty(found),'JMPMD:CJKFont','Install a Chinese font, or pass Font explicitly. Bundled PDFs remain viewable without MATLAB fonts.');
    fontName=choices{found};
    % Windows listfonts may return the localized family name for SimSun.
    if strcmp(fontName,'宋体'),fontName='SimSun';end
end
if ~exist(outDir,'dir'),mkdir(outDir);end
outDir=char(java.io.File(outDir).getCanonicalPath());
published=char(java.io.File(fullfile(root,'results','published')).getCanonicalPath());
assert(~strcmpi(outDir,published)&&~startsWith(lower(outDir),[lower(published) filesep]), ...
    'JMPMD:PublishedReadOnly','Export outside results/published.');
tables=struct;
for name={'summary','convergence'}
    file=fullfile(dataDir,[name{1} '.csv']);
    if exist(file,'file'),tables.(name{1})=readtable(file,'TextType','string','VariableNamingRule','preserve');end
end
spec=jsondecode(fileread(fullfile(root,'metadata','figure_layout.json')));ready=cell(0,1);skipped={};
for fi=1:numel(spec)
    f=spec(fi);reason='';
    for pi=1:numel(f.panels)
        for j=1:numel(f.panels(pi).curves)
            c=f.panels(pi).curves(j);
            if ~isfield(tables,c.table),reason=['Missing ' c.table '.csv'];break;end
            tab=tables.(c.table);mask=tab.Kind==string(c.kind);
            if strcmp(c.table,'summary')
                mask=mask&tab.Group==string(c.group)&tab.Series==string(c.series);xfield='X';
            else
                mask=mask&tab.SNR==c.snr;xfield='Iteration';
            end
            rows=tab(mask,:);x=c.expected_x(:);y=nan(size(x));err=[];
            if strcmp(c.table,'convergence') && ismember('Profile',rows.Properties.VariableNames)&&~isempty(rows)&&all(rows.Profile=="smoke")
                x=x(x<=max(rows.Iteration));
            end
            for k=1:numel(x)
                idx=find(rows.(xfield)==x(k));
                if numel(idx)~=1,reason=sprintf('Missing or duplicate %s/%s/%s x=%g',c.group,c.kind,c.series,x(k));break;end
                if ismember('Pending',rows.Properties.VariableNames)&&rows.Pending(idx)>0,reason='Requested trials are incomplete';break;end
                y(k)=rows.(c.metric)(idx);
                if ~isempty(c.sem_metric),err(k,1)=rows.(c.sem_metric)(idx);end %#ok<AGROW>
            end
            y=y(1:numel(x));
            if isempty(x)||any(~isfinite(y))||any(~isfinite(err))||(f.panels(pi).log&&any(y<=0))
                reason='Nonfinite/nonpositive log metric (retained in CSV; no clipping or pseudo-count)';
            end
            if ~isempty(reason),break;end
            f.panels(pi).curves(j).x=x;f.panels(pi).curves(j).y=y;f.panels(pi).curves(j).error=err;
        end
        if ~isempty(reason),break;end
    end
    if isempty(reason),ready{end+1}=f;else,skipped{end+1}=struct('figure',f.name,'reason',reason);end %#ok<AGROW>
end
if isempty(ready)
    audit=struct('figures',0,'curves',0,'points',0);
else
    audit=render_figures([ready{:}],outDir,fontName);
end
audit.skipped=skipped;audit.source_data=char(dataDir);audit.font=fontName;
atomic_json(fullfile(outDir,'plot_audit.json'),audit);
fprintf('Rendered %d main figures; skipped %d. See plot_audit.json.\n',audit.figures,numel(skipped));
end
