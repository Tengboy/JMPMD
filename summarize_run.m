function [summary,convergence]=summarize_run(out)
%SUMMARIZE_RUN Rebuild lightweight CSV from committed batches only.
% Means/SEMs retain failed samples, NaN, Inf and extreme errors. No filtering.
setup_jmpmd;d=load(fullfile(out,'manifest.mat'));
[done,stats]=scan_chunks(out,d.conditions,d.signature);
summary=table;convergence=table;
for ci=1:numel(d.conditions)
    c=d.conditions(ci);data=nan(c.n_trials,43);available=false(c.n_trials,1);
    nTrace=min(c.trace_trials,done(ci));T=c.cfg.max_iterations;
    if nTrace>0,traces=nan(T,5,nTrace);end
    files=dir(fullfile(out,'chunks',sprintf('c%03d_*.mat',ci)));
    for k=1:numel(files)
        b=load(fullfile(files(k).folder,files(k).name));data(b.trial_ids,:)=b.records;available(b.trial_ids)=true;
        for j=1:numel(b.trial_ids)
            id=b.trial_ids(j);
            if id<=nTrace&&~isempty(b.traces{j})
                assert(isequal(size(b.traces{j}),[T 5]),'JMPMD:TraceShape','Invalid convergence trace.');
                traces(:,:,id)=b.traces{j};
            end
        end
    end
    v=d.views([d.views.condition_id]==ci);[v.condition_id]=deal(1);one=c;one.id=1;
    local=struct('conditions',one,'views',v,'records',{{data}},'done',{{available}});
    rows=summarize_records(local);rows.Condition(:)=ci;
    rows.SourceRoot(:)="local_run";rows.SourceCondition(:)=ci;rows.SourceSignature(:)=string(d.signature);
    rows.InitMode=repmat(string(c.cfg.init_mode),height(rows),1);rows.AngleSolver=repmat("local",height(rows),1);
    rows.GainPilotInit=true(height(rows),1);rows.Profile=repmat(string(d.profile),height(rows),1);
    summary=[summary;rows]; %#ok<AGROW>
    if nTrace>0
        mu=mean(traces,3);se=std(traces,0,3)/sqrt(nTrace);
        if nTrace<2,se(:)=NaN;end
        active=sum(data(1:nTrace,8).'>=(1:T)',2);
        Kind=repmat(string(c.cfg.kind),T,1);SNR=repmat(c.cfg.snr_db,T,1);Iteration=(1:T)';N=repmat(nTrace,T,1);
        tab=table(Kind,SNR,Iteration,N,mu(:,1),mu(:,2),mu(:,3),mu(:,4),mu(:,5), ...
            se(:,2),se(:,3),se(:,5),active,'VariableNames',{'Kind','SNR','Iteration','N','G_MSE','G_NMSE', ...
            'X_MSE','Theta','Position','G_NMSE_SEM','X_MSE_SEM','Position_SEM','Active'});
        tab.SourceCondition=repmat(ci,T,1);tab.Profile=repmat(string(d.profile),T,1);
        tab.Pending=repmat(c.trace_trials-nTrace,T,1);convergence=[convergence;tab]; %#ok<AGROW>
    end
end
writetable(summary,fullfile(out,'summary.csv'),'Encoding','UTF-8');
if ~isempty(convergence),writetable(convergence,fullfile(out,'convergence.csv'),'Encoding','UTF-8');end
atomic_json(fullfile(out,'summary_audit.json'),struct('source_signature',d.signature,'profile',d.profile, ...
    'completed_trials',sum(done),'target_trials',sum([d.conditions.n_trials]),'batches',stats.batches, ...
    'summary_rows',height(summary),'convergence_rows',height(convergence), ...
    'failed_trials',stats.failed,'reference_failed',stats.reference_failed, ...
    'mean_policy','All committed requested trials, including NaN/Inf/failures; no tail removal.'));
end
