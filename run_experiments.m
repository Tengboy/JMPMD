function status=run_experiments(varargin)
%RUN_EXPERIMENTS Run or resume the selected main-figure Monte Carlo trials.
% run_experiments('Profile','smoke')                       % 440 short trials
% run_experiments('Profile','full','Workers',12)            % 1,300,000 trials
% run_experiments('Profile','full','Groups',{'B3'},'Kinds',{'Q'})
% No automatic/background scheduler is installed by this function.
root=setup_jmpmd;
p=inputParser;p.addParameter('Profile','smoke',@(x)ismember(string(x),["smoke","full"]));
p.addParameter('Groups',{'A1','A2','A3','A4','B1','B3'});p.addParameter('Kinds',{'Q','G'});
p.addParameter('Workers',0,@(x)isscalar(x)&&isfinite(x)&&x>=0&&x==fix(x));
p.addParameter('BatchSize',250,@(x)isscalar(x)&&isfinite(x)&&x>=1&&x==fix(x));
p.addParameter('MaxBatches',Inf,@(x)isscalar(x)&&x>=1&&x==fix(x));p.addParameter('Output','');p.parse(varargin{:});opt=p.Results;
[conditions,views]=experiment_catalog(opt.Groups,opt.Kinds,opt.Profile);columns=record_columns;
out=char(opt.Output);
if isempty(out),out=fullfile(root,'runs',char(opt.Profile));end
if ~exist(out,'dir'),mkdir(out);end
out=char(java.io.File(out).getCanonicalPath());
published=char(java.io.File(fullfile(root,'results','published')).getCanonicalPath());
assert(~strcmpi(out,published)&&~startsWith(lower(out),[lower(published) filesep]), ...
    'JMPMD:PublishedReadOnly','Choose an output outside results/published.');
if ~exist(fullfile(out,'chunks'),'dir'),mkdir(fullfile(out,'chunks'));end
lockPath=fullfile(out,'RUNNING.lock');lock=java.io.File(lockPath);
assert(lock.createNewFile(),'JMPMD:RunLocked','RUNNING.lock exists. Verify its process before removing a stale lock; see docs/REPRODUCING.md.');
cleanup=onCleanup(@()remove_lock(lockPath)); %#ok<NASGU>
pid=feature('getpid');fid=fopen(lockPath,'w');fprintf(fid,'pid=%d\nhost=%s\n',pid,char(java.net.InetAddress.getLocalHost.getHostName));fclose(fid);
[signature,sourceManifest]=source_signature(root);manifestPath=fullfile(out,'manifest.mat');
if exist(manifestPath,'file')
    previous=load(manifestPath);
    assert(strcmp(previous.signature,signature),'JMPMD:ChangedSource','Source changed; use a new output directory.');
    assert(isequaln(previous.conditions,conditions)&&isequaln(previous.views,views), ...
        'JMPMD:ChangedConfig','Selection/profile changed; use a new output directory.');
    assert(strcmp(previous.matlab_version,version),'JMPMD:ChangedMATLAB','MATLAB version changed; use a new output directory.');
    started_at=previous.started_at;
else
    matlab_version=version;started_at=datestr(now,31);profile=char(opt.Profile);
    save([manifestPath '.pending'],'conditions','views','columns','signature','sourceManifest', ...
        'matlab_version','started_at','profile','-v7');
    [ok,msg]=movefile([manifestPath '.pending'],manifestPath);assert(ok,msg);
    atomic_json(fullfile(out,'manifest.json'),struct('conditions',conditions,'views',views, ...
        'columns',{columns},'signature',signature,'matlab_version',matlab_version,'profile',profile,'seed_base',20260927));
end
atomic_json(fullfile(out,'process.json'),struct('pid',pid,'host',char(java.net.InetAddress.getLocalHost.getHostName), ...
    'started_at',datestr(now,31),'profile',opt.Profile,'signature',signature));
[done,stats]=scan_chunks(out,conditions,signature);write_status('starting');
threadCount=maxNumCompThreads;maxNumCompThreads(1);threadCleanup=onCleanup(@()maxNumCompThreads(threadCount)); %#ok<NASGU>
pool=[];ownPool=false;
if opt.Workers>0 && any(done<[conditions.n_trials])
    assert(license('test','Distrib_Computing_Toolbox'),'JMPMD:ParallelToolbox','Use Workers=0 for base-MATLAB serial execution.');
    pool=gcp('nocreate');
    if isempty(pool)
        cluster=parcluster('local');storage=fullfile(out,'parallel_jobs');if ~exist(storage,'dir'),mkdir(storage);end
        cluster.JobStorageLocation=storage;pool=parpool(cluster,opt.Workers);ownPool=true;
    else
        assert(pool.NumWorkers==opt.Workers,'JMPMD:ExistingPool','An existing pool has a different size; choose Workers=0 or its actual size.');
    end
    if ownPool,poolCleanup=onCleanup(@()delete(pool));end %#ok<NASGU>
    f=parfevalOnAll(pool,@maxNumCompThreads,1,1);wait(f);
end
batches=0;
for ci=1:numel(conditions)
    condition=conditions(ci);
    while done(ci)<condition.n_trials
        trial_ids=done(ci)+(1:min(opt.BatchSize,condition.n_trials-done(ci)));n=numel(trial_ids);
        records=nan(n,43);traces=cell(n,1);errors=cell(n,1);reference_errors=cell(n,1);warning_ids=cell(n,1);timer=tic;
        if isempty(pool)
            for j=1:n,[records(j,:),traces{j},errors{j},reference_errors{j},warning_ids{j}]=simulate_trial(condition,trial_ids(j));end
        else
            parfor j=1:n,[records(j,:),traces{j},errors{j},reference_errors{j},warning_ids{j}]=simulate_trial(condition,trial_ids(j));end
        end
        wall_seconds=toc(timer);condition_id=ci;seeds=20260927+trial_ids;
        file=fullfile(out,'chunks',sprintf('c%03d_%06d_%06d.mat',ci,trial_ids(1),trial_ids(end)));
        assert(~exist(file,'file'),'JMPMD:DuplicateChunk','Will not overwrite committed trials.');
        save([file '.pending'],'condition_id','trial_ids','seeds','records','traces','errors', ...
            'reference_errors','warning_ids','wall_seconds','signature','-v7');
        [ok,msg]=movefile([file '.pending'],file);assert(ok,msg);done(ci)=trial_ids(end);
        stats.batches=stats.batches+1;stats.wall_seconds=stats.wall_seconds+wall_seconds;
        stats.failed=stats.failed+nnz(records(:,21)==0);stats.reference_failed=stats.reference_failed+nnz(records(:,42)==0);
        batches=batches+1;write_status('running');
        fprintf('%s %s c%03d %c/%c %d/%d; total %d/%d; failed %d\n',datestr(now,31), ...
            opt.Profile,ci,condition.cfg.kind,condition.cfg.method,done(ci),condition.n_trials,sum(done),sum([conditions.n_trials]),stats.failed);
        if (batches>=opt.MaxBatches||exist(fullfile(out,'STOP_AFTER_CHECKPOINT'),'file'))&&any(done<[conditions.n_trials])
            summarize_run(out);write_status('stopped_after_checkpoint');return;
        end
    end
end
summarize_run(out);stage='completed';if stats.failed||stats.reference_failed,stage='completed_with_failures';end
write_status(stage);
fprintf('Completed. CSV results: %s\n',out);
    function write_status(stage)
        status=struct('stage',stage,'profile',opt.Profile,'started_at',started_at,'updated_at',datestr(now,31), ...
            'pid',pid,'signature',signature,'completed',sum(done),'target',sum([conditions.n_trials]), ...
            'conditions',numel(conditions),'views',numel(views),'completed_by_condition',done, ...
            'committed_batches',stats.batches,'batch_wall_seconds',stats.wall_seconds, ...
            'failed_trials',stats.failed,'reference_failed',stats.reference_failed);
        atomic_json(fullfile(out,'status.json'),status);
    end
end
function remove_lock(p),if exist(p,'file'),delete(p);end,end
