function audit=verify_run(out)
%VERIFY_RUN Validate checkpoints and summary availability (no re-simulation).
root=setup_jmpmd;d=load(fullfile(out,'manifest.mat'));[signature,~]=source_signature(root);
assert(strcmp(signature,d.signature),'JMPMD:ChangedSource','Source differs from run manifest.');
[done,stats]=scan_chunks(out,d.conditions,d.signature);
[summary,conv]=summarize_run(out);
assert(all(done==[d.conditions.n_trials]),'JMPMD:Incomplete','Not all trials are committed.');
assert(height(summary)==numel(d.views)&&all(summary.Pending==0));
cfgs=[d.conditions.cfg];assert(all([cfgs.strict]),'JMPMD:StrictOnly','Non-strict config.');
traceTarget=sum([d.conditions.trace_trials]);
if traceTarget>0
    traceCfg=d.conditions([d.conditions.trace_trials]>0);
    assert(height(conv)==sum(arrayfun(@(x)x.cfg.max_iterations,traceCfg)));
end
audit=struct('structure_valid',true,'profile',d.profile,'conditions',numel(d.conditions),'views',numel(d.views), ...
    'trials',sum(done),'batches',stats.batches,'trajectory_trials',traceTarget, ...
    'failed_trials',stats.failed,'reference_failed',stats.reference_failed, ...
    'primary_nonfinite_across_views',sum(summary.PrimaryNonfinite), ...
    'note','View totals include shared conditions; numerical failures are preserved and do not erase structure validation.');
atomic_json(fullfile(out,'verification.json'),audit);
end
