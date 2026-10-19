function summary=summarize_records(d)
% Ordinary means include every requested result. No omitnan or tail filtering.
rows=cell(1,numel(d.views));
for iv=1:numel(d.views)
    v=d.views(iv);c=d.conditions(v.condition_id);cfg=c.cfg;
    selected=1:v.n_use;available=d.done{v.condition_id}(selected);data=d.records{v.condition_id}(selected(available),:);
    r=struct('Group',string(v.group),'Condition',v.condition_id,'Kind',string(cfg.kind), ...
        'Method',string(cfg.method),'Series',string(v.series),'X',v.x,'SNR',cfg.snr_db, ...
        'NoiseMode',string(cfg.noise_mode),'NB',cfg.BS_size,'M',cfg.RIS_m,'NR',cfg.RIS_size_x*cfg.RIS_size_y, ...
        'Tp',cfg.T_p,'Td',cfg.T-cfg.T_p,'NLoS',string(cfg.nlos_model),'Kappa',cfg.kappa_db,'Strict',cfg.strict, ...
        'Target',v.n_use,'N',size(data,1),'Pending',v.n_use-size(data,1));
    r.Failed=nnz(data(:,21)==0);r.AtCap=nnz(data(:,10)==1);r.Unconverged=nnz(data(:,9)==0);
    r.TaylorCalls=sum_known(data(:,12));r.TaylorFailedCalls=sum_known(data(:,13));r.TaylorAffectedTrials=nnz(data(:,13)>0);
    r.ReferenceFailed=nnz(data(:,42)==0);r.ReferenceMissing=0;
    if c.need_reference,r.ReferenceMissing=nnz(~isfinite(data(:,42)));end
    r.LargePositionTrials=nnz(any(data(:,5:7)>1,2));r.LargePositionUsers=nnz(data(:,5:7)>1);
    r.UnknownPositionUsers=nnz(~isfinite(data(:,5:7)));r.PrimaryNonfinite=nnz(~isfinite(data(:,1:3)));
    r.Theta=average(data,1);r.Position=average(data,2);r.Data=average(data,3);r.ThetaPost=average(data,4);
    r.ThetaSEM=sem(data,1);r.PositionSEM=sem(data,2);r.DataSEM=sem(data,3);
    r.Iterations=average(data,8);r.EstimatorSeconds=average(data,19);r.TotalSeconds=average(data,20);
    r.ActualSNR=average(data,22);r.InputKappa=average(data,23);r.ProjectedKappa=average(data,24);
    r.PositionBound=average(data,33);r.ThetaBound=average(data,34);r.DataBound=average(data,35);
    r.CSISER=average(data,36);r.CSIFullSER=average(data,37);r.BoundSeconds=average(data,38);
    r.MaxBoundResidual=max_or_nan(data(:,41));r.MaxWorkerMemoryGB=max_or_nan(data(:,43))/2^30;
    r.SymbolErrors=NaN;r.SymbolDecisions=NaN;
    if cfg.kind=='Q',r.SymbolErrors=round(sum(data(:,3))*cfg.User_k*(cfg.T-cfg.T_p));r.SymbolDecisions=r.N*cfg.User_k*(cfg.T-cfg.T_p);end
    r.SourceRoot="";r.SourceCondition=NaN;r.SourceSignature="";
    if isfield(d,'provenance')&&~isempty(d.provenance)
        p=d.provenance(v.condition_id);r.SourceRoot=p.SourceRoot;r.SourceCondition=p.SourceCondition;r.SourceSignature=p.SourceSignature;
    end
    if r.Pending>0,r.Theta=NaN;r.Position=NaN;r.Data=NaN;end
    rows{iv}=r;
end
summary=struct2table([rows{:}]);
end
function value=average(data,col),if isempty(data),value=NaN;else,value=mean(data(:,col));end,end
function value=sem(data,col),if size(data,1)<2,value=NaN;else,value=std(data(:,col))/sqrt(size(data,1));end,end
function value=sum_known(x),value=sum(x(isfinite(x)));end
function value=max_or_nan(x),if isempty(x),value=NaN;else,value=max(x);end,end
