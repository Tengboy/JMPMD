function [conditions,views]=experiment_catalog(groups,kinds,profile)
%EXPERIMENT_CATALOG Exact current main-figure settings, deduplicated by cfg.
% J=JMPMD; P=pilot-only; W=no position feedback. Q=QPSK, G=Gaussian.
if nargin<1,groups={'A1','A2','A3','A4','B1','B3'};end
if nargin<2,kinds={'Q','G'};end
if nargin<3,profile='full';end
groups=cellstr(upper(string(groups)));kinds=cellstr(upper(string(kinds)));
assert(all(ismember(groups,{'A1','A2','A3','A4','B1','B3'}))&&~isempty(groups),'JMPMD:Groups','Unknown/empty group selection.');
assert(all(ismember(kinds,{'Q','G'}))&&~isempty(kinds),'JMPMD:Kinds','Kinds must be Q and/or G.');
assert(ismember(profile,{'full','smoke'}),'JMPMD:Profile','Choose full or smoke.');
conditions=struct('id',{},'key',{},'cfg',{},'n_trials',{},'need_reference',{},'trace_trials',{});
views=struct('group',{},'condition_id',{},'n_use',{},'x',{},'series',{});
index=containers.Map('KeyType','char','ValueType','double');
for kind='QG'
    if ~ismember(string(kind),string(kinds)),continue;end
    for snr=-5:10
        c=defaults(kind);c.snr_db=snr;n=5000;
        if kind=='Q'&&snr>=7,n=15000;end
        for method='JPW',c.method=method;register('A1',c,n,snr,method,method=='J',0);end
        c.method='J';c.nlos_model='mixed20';
        for ka=[10 15],c.kappa_db=ka;register('A2',c,n,snr,sprintf('K%d',ka),false,0);end
    end
    for nb=16:4:40,c=defaults(kind);c.BS_size=nb;register('A3',c,5000,nb,'J',true,0);end
    for snr=[0 6],c=defaults(kind);c.snr_db=snr;register('A4',c,5000,snr,'J',false,5000);end
    for tp=[5 10 15 20]
        for method='JPW',c=defaults(kind);c.method=method;c.T_p=tp;register('B1',c,5000,tp,method,false,0);end
    end
    for side=8:2:16
        for method='JPW'
            c=defaults(kind);c.method=method;c.noise_mode='fixed';
            c.RIS_size_x=side;c.RIS_size_y=side;c.init_mode='pilot_local';
            register('B3',c,5000,side,method,false,0);
        end
    end
end
    function register(group,cfg,n,x,series,reference,traceN)
        if ~ismember(group,groups),return;end
        if strcmp(profile,'smoke')
            n=2;cfg.max_iterations=4;
            if traceN>0,traceN=2;end
        end
        key=jsonencode(cfg);
        if isKey(index,key)
            id=index(key);conditions(id).n_trials=max(conditions(id).n_trials,n);
            conditions(id).need_reference=conditions(id).need_reference||reference;
            conditions(id).trace_trials=max(conditions(id).trace_trials,traceN);
        else
            id=numel(conditions)+1;index(key)=id;
            conditions(id)=struct('id',id,'key',key,'cfg',cfg,'n_trials',n,'need_reference',reference,'trace_trials',traceN);
        end
        views(end+1)=struct('group',group,'condition_id',id,'n_use',n,'x',x,'series',series);
    end
end

function cfg=defaults(kind)
cfg=struct('kind',kind,'method','J','strict',true,'snr_db',6, ...
    'noise_mode','snr','noise_dbm',-125,'Tx_power_dB',25,'f_c',28e9, ...
    'BS_size',32,'RIS_m',3,'RIS_size_x',10,'RIS_size_y',10,'User_k',3, ...
    'T',40,'T_p',10,'pos_var',0.1,'max_iterations',200, ...
    'stop_tolerance',1e-5,'stop_stable_iterations',2, ...
    'nlos_model','los','kappa_db',Inf,'init_mode','prior', ...
    'gain_pilot_init',true,'angle_solver','local');
end
