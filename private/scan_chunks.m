function [done,stats]=scan_chunks(out,conditions,signature)
% Validate every committed batch before resume, never repair or skip it.
done=zeros(1,numel(conditions));stats=struct('batches',0,'wall_seconds',0,'failed',0,'reference_failed',0);
paths=dir(fullfile(out,'chunks','c*.mat'));
for k=1:numel(paths)
    b=load(fullfile(paths(k).folder,paths(k).name));ci=b.condition_id;
    assert(ci>=1&&ci<=numel(conditions)&&ci==fix(ci),'JMPMD:Condition','Invalid chunk condition.');
    assert(strcmp(b.signature,signature),'JMPMD:Signature','Chunk source signature differs.');
    assert(~isempty(b.trial_ids)&&isequal(b.trial_ids,done(ci)+(1:numel(b.trial_ids))), ...
        'JMPMD:ChunkGap','Duplicate, missing or out-of-order committed trials.');
    assert(isequal(b.seeds,20260927+b.trial_ids)&&b.trial_ids(end)<=conditions(ci).n_trials, ...
        'JMPMD:Seeds','Seed/trial range mismatch.');
    assert(isequal(size(b.records),[numel(b.trial_ids) 43]),'JMPMD:Records','Invalid record shape.');
    assert(numel(b.traces)==numel(b.trial_ids),'JMPMD:Traces','Invalid trace count.');
    done(ci)=b.trial_ids(end);stats.batches=stats.batches+1;
    stats.wall_seconds=stats.wall_seconds+b.wall_seconds;
    stats.failed=stats.failed+nnz(b.records(:,21)==0);
    stats.reference_failed=stats.reference_failed+nnz(b.records(:,42)==0);
end
end
