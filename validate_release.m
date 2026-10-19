function report=validate_release(varargin)
%VALIDATE_RELEASE Check math, catalog, golden samples, checkpoint resume.
% Default takes a few minutes, never runs a full Monte Carlo experiment.
% validate_release('Output',fullfile(tempdir,'jmpmd_validation'))
root=setup_jmpmd;p=inputParser;p.addParameter('Output',fullfile(root,'validation_run'));
p.addParameter('Golden',true,@islogical);p.parse(varargin{:});out=char(p.Results.Output);
if ~exist(out,'dir'),mkdir(out);end
assert(~exist(fullfile(out,'resume','manifest.mat'),'file'),'JMPMD:FreshValidation','Use a fresh validation output directory.');
addpath(fullfile(root,'tests'));test_strict_moments(fullfile(root,'src'));
[c,v]=experiment_catalog;cfgs=[c.cfg];
assert(numel(c)==220&&numel(v)==232&&sum([c.n_trials])==1300000);
assert(sum([c.trace_trials])==20000&&sum([c.need_reference])==44);
assert(sum([c(strcmp({cfgs.kind},'Q')).n_trials])==750000);
assert(sum([c(strcmp({cfgs.kind},'G')).n_trials])==550000);
assert(all([cfgs.strict])&&all(strcmp({cfgs.angle_solver},'local')));
% Independently establish physical common observations and SNR/kappa.
a=c(1).cfg;s=make_scene(a,20260928);a.method='P';sp=make_scene(a,20260928);
assert(isequaln(s.Y,sp.Y));a.nlos_model='mixed20';a.kappa_db=10;sn=make_scene(a,20260928);
assert(isequaln(s.Y_los,sn.Y_los)&&isequaln(s.noise,sn.noise)&&isequaln(s.X,sn.X));
assert(abs(sn.actual_snr_db-a.snr_db)<1e-12&&abs(sn.nlos.input_kappa_db-10)<1e-12);
% A failed sample must poison ordinary means, not disappear.
one=c(1);one.n_trials=2;one.trace_trials=0;one.need_reference=false;
data=ones(2,43);data(2,1:3)=NaN;data(2,21)=0;
vv=struct('group','A1','condition_id',1,'n_use',2,'x',-5,'series','J');
d=struct('conditions',one,'views',vv,'records',{{data}},'done',{{true(2,1)}});
t=summarize_records(d);assert(t.N==2&&t.Failed==1&&t.PrimaryNonfinite==3&&isnan(t.Theta)&&isnan(t.Position)&&isnan(t.Data));
goldenTrials=0;maxScaledError=0;
if p.Results.Golden
    g=readtable(fullfile(root,'tests','golden_samples.csv'),'TextType','string');
    names=record_columns;cols=[1:18 21:37 39:42];
    for k=1:height(g)
        condition=c(g.ReleaseCondition(k));[actual,~,err,referr]=simulate_trial(condition,g.Trial(k));
        assert(isempty(err)&&isempty(referr),'JMPMD:GoldenExecution','Golden sample failed.');
        expected=table2array(g(k,names(cols)));
        assert(isequal(isnan(actual(cols)),isnan(expected))&&isequal(isinf(actual(cols)),isinf(expected)));
        mask=isfinite(expected);a=actual(cols);delta=abs(a(mask)-expected(mask))./max(abs(expected(mask)),1e-12);
        maxScaledError=max(maxScaledError,max(delta));
        assert(all(delta<2e-8),'JMPMD:GoldenMismatch','Frozen sample differs (row %d).',k);
        goldenTrials=goldenTrials+1;
    end
end
% Full-length A4 traces use held final values after the estimator stops.
[traceConditions,~]=experiment_catalog({'A4'},{'Q','G'},'full');
for k=1:numel(traceConditions)
    [r,tr,err,referr]=simulate_trial(traceConditions(k),1);
    assert(isempty(err)&&isempty(referr)&&r(21)==1&&isequal(size(tr),[200 5]));
    assert(abs(tr(end,5)-r(2))<1e-12);
    finalIter=r(8);assert(finalIter>=1&&finalIter<=200);
    tail=repmat(tr(finalIter,:),201-finalIter,1);
    assert(isequaln(tr(finalIter:end,:),tail),'JMPMD:TraceExtension','Convergence trace is not held at its final value.');
end
% Stop after one committed trial, then resume. Changing cfg must be rejected.
resume=fullfile(out,'resume');args={'Profile','smoke','Groups',{'A4'},'Workers',0,'BatchSize',1,'Output',resume};
first=run_experiments(args{:},'MaxBatches',1);assert(first.completed==1);
file=dir(fullfile(resume,'chunks','*.mat'));fid=fopen(fullfile(file.folder,file.name),'rb');before=fread(fid,Inf,'*uint8');fclose(fid);
try
    run_experiments('Profile','smoke','Groups',{'A3'},'Output',resume);error('Expected config rejection.');
catch e
    assert(strcmp(e.identifier,'JMPMD:ChangedConfig'));
end
status=run_experiments(args{:});assert(status.completed==8&&status.failed_trials==0);
fid=fopen(fullfile(file.folder,file.name),'rb');after=fread(fid,Inf,'*uint8');fclose(fid);assert(isequal(before,after));
again=run_experiments(args{:});assert(again.committed_batches==status.committed_batches);
verify=verify_run(resume);assert(verify.trials==8&&verify.trajectory_trials==8);
% Independently inject a duplicate only into a copied validation fixture.
duplicate=fullfile(out,'duplicate_fixture');copyfile(resume,duplicate);
copyfile(fullfile(file.folder,file.name),fullfile(duplicate,'chunks','c001_000001_duplicate.mat'));
try
    verify_run(duplicate);error('Expected duplicate rejection.');
catch e
    assert(strcmp(e.identifier,'JMPMD:ChunkGap'));
end
% Simulate an existing lock: a second launcher must not enter that run.
fid=fopen(fullfile(resume,'RUNNING.lock'),'w');fprintf(fid,'validation fixture');fclose(fid);
try
    run_experiments(args{:});error('Expected lock rejection.');
catch e
    delete(fullfile(resume,'RUNNING.lock'));assert(strcmp(e.identifier,'JMPMD:RunLocked'));
end
report=struct('passed',true,'catalog_conditions',220,'catalog_views',232,'full_trials',1300000, ...
    'golden_trials',goldenTrials,'max_scaled_golden_error',maxScaledError, ...
    'strict_moment_checks',true,'scene_checks',true,'failed_mean_retained',true, ...
    'full_A4_trace_checks',4,'resume_committed_unchanged',true,'changed_config_rejected',true, ...
    'duplicate_chunk_rejected',true,'existing_lock_rejected',true,'matlab_version',version);
atomic_json(fullfile(out,'validation.json'),report);disp(report);
end
