function [values,trace,errorText,referenceError,warningId]=simulate_trial(condition,trial)
%SIMULATE_TRIAL A self-contained trial; no historical MAT caches are read.
% All methods use seed=20260927+trial. References are computed, not imported.
started=tic;values=nan(1,43);trace=[];errorText='';referenceError='';warningId='';
cfg=condition.cfg;seed=20260927+trial;lastwarn('');
assert(cfg.strict&&strcmp(cfg.angle_solver,'local'),'JMPMD:StrictOnly','Release catalog uses strict moments and local solver.');
try
    scene=make_scene(cfg,seed);
    values(22:32)=[scene.actual_snr_db,scene.nlos.input_kappa_db,scene.nlos.projected_kappa_db, ...
        mean(abs(scene.Y_los(:)).^2),scene.interference_power,scene.noise_power, ...
        scene.nlos.mixed_cascade_power,scene.nlos.mixed_direct_power,scene.nlos.cross_power, ...
        sum(real(scene.Y(:))),sum(imag(scene.Y(:)))];
    iniSeconds=0;
    if strcmp(cfg.init_mode,'pilot_local')
        initTimer=tic;
        names={'X_p','Y','A_BS','Omega','effective_noise_power','User_pos_c','RIS_pos','RIS_dir_x','RIS_dir_y'};
        observed=struct;for k=1:numel(names),observed.(names{k})=scene.(names{k});end
        cfg.initial_theta=pilot_initial_angles(cfg,observed);iniSeconds=toc(initTimer);
    end
    estimator=str2func(['estimate_' lower(cfg.kind) '_' lower(cfg.method)]);
    if cfg.kind=='G'
        if cfg.method=='J',estimator=@estimate_g_j_warm;end
        if cfg.method=='W',estimator=@estimate_g_w_warm;end
    end
    result=estimator(cfg,scene,trial<=condition.trace_trials);
    values(1:18)=[result.theta_mse,result.position_mse,result.data_error,result.theta_post_mse, ...
        result.position_user_error,result.iterations,result.converged,result.cap,result.last_change, ...
        result.localization_calls,result.localization_failures,result.localization_reasons];
    values(19)=result.elapsed+iniSeconds;values(21)=result.finite;trace=result.trace;
catch problem
    values(21)=0;errorText=getReport(problem,'extended','hyperlinks','off');
end
% Keep reference failure independent of estimator failure whenever a scene exists.
if condition.need_reference
    referenceTimer=tic;
    try
        ref=suite_reference(cfg,scene);
        values(33:42)=[ref.position_bound,ref.theta_bound,ref.data_bound,ref.csi_ser,ref.csi_full_ser, ...
            toc(referenceTimer),ref.fim_rcond,ref.scaled_fim_rcond,ref.inverse_residual,ref.finite];
    catch problem
        referenceError=getReport(problem,'extended','hyperlinks','off');values(38)=toc(referenceTimer);values(42)=0;
    end
end
[~,warningId]=lastwarn;values(20)=toc(started);
try,s=memory;values(43)=s.MemUsedMATLAB;catch,end
end
