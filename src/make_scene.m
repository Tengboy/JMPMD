function s=make_scene(cfg,seed)
% One observation shared by J/P/W. Separate deterministic random streams
% keep LoS/symbols/probes/noise unchanged when NLoS or estimator changes.
s=cfg; s.seed=seed; s.lambda=3e8/cfg.f_c;
s.Tx_power=0.001*10^(cfg.Tx_power_dB/10);
s.T_d=cfg.T-cfg.T_p; s.RIS_size=cfg.RIS_size_x*cfg.RIS_size_y;
s.BS_pos=[-20;0;0]; s.BS_dir=[0;1;0];
s.RIS_pos=[0 0 20;20 -20 0;10 10 5]; s.RIS_pos=s.RIS_pos(:,1:cfg.RIS_m);
s.RIS_dir_x=[1 -1 0;0 0 1;0 0 0]; s.RIS_dir_x=s.RIS_dir_x(:,1:cfg.RIS_m);
s.RIS_dir_y=repmat([0;0;1],1,cfg.RIS_m);
s.User_pos_c=[-5 10 10;0 10 -10;3.5 1 1];
s.RIS_para=struct('RIS_dir_x',s.RIS_dir_x,'RIS_dir_y',s.RIS_dir_y, ...
    'RIS_size_x',cfg.RIS_size_x,'RIS_size_y',cfg.RIS_size_y,'RIS_pos',s.RIS_pos);
r=RandStream('mt19937ar','Seed',seed);
s.User_pos=s.User_pos_c+sqrt(cfg.pos_var)*randn(r,3,cfg.User_k);
s.d_B=vecnorm(s.BS_pos-s.RIS_pos).';
s.phi_RIS_x=pi*diag(s.RIS_dir_x.'*(s.BS_pos-s.RIS_pos))./s.d_B;
s.phi_RIS_y=pi*diag(s.RIS_dir_y.'*(s.BS_pos-s.RIS_pos))./s.d_B;
s.th_BS=pi*(s.BS_dir.'*(s.RIS_pos-s.BS_pos)).'./s.d_B;
s.A_BS=array(s.th_BS,cfg.BS_size);
s.th_RIS_x=zeros(cfg.User_k,cfg.RIS_m); s.th_RIS_y=s.th_RIS_x;
s.alpha=zeros(cfg.User_k,cfg.RIS_m); s.d_U=s.alpha;
for m=1:cfg.RIS_m
    offset=s.User_pos-s.RIS_pos(:,m); distance=vecnorm(offset);
    s.d_U(:,m)=distance.';
    s.th_RIS_x(:,m)=(pi*s.RIS_dir_x(:,m).'*offset./distance).';
    s.th_RIS_y(:,m)=(pi*s.RIS_dir_y(:,m).'*offset./distance).';
    s.alpha(:,m)=sqrt(s.Tx_power)*s.lambda^2*exp(-1j*2*pi*(s.d_B(m)+distance.')/s.lambda) ...
        ./(16*pi^2*s.d_B(m)*distance.');
end
s.th_RIS=cat(3,s.th_RIS_x,s.th_RIS_y);
s.constellation=[1+1j -1+1j -1-1j 1-1j]/sqrt(2);
r=RandStream('mt19937ar','Seed',seed+2000000);
if cfg.kind=='Q'
    symbols=randi(r,4,cfg.T,cfg.User_k); s.X=s.constellation(symbols);
    s.data_ind_p=symbols(1:cfg.T_p,:); s.data_ind_d=symbols(cfg.T_p+1:end,:);
else
    s.X=(randn(r,cfg.T,cfg.User_k)+1j*randn(r,cfg.T,cfg.User_k))/sqrt(2);
    s.data_ind_p=[]; s.data_ind_d=[];
end
s.X_p=s.X(1:cfg.T_p,:); s.X_d=s.X(cfg.T_p+1:end,:);
r=RandStream('mt19937ar','Seed',seed+4000000);
s.Omega=exp(1j*2*pi*rand(r,s.RIS_size,cfg.T,cfg.RIS_m));
s.rho=zeros(cfg.T,cfg.RIS_m); s.H=zeros(cfg.RIS_m,cfg.User_k,s.T_d);
for m=1:cfg.RIS_m
    s.Omega(:,:,m)=s.Omega(:,:,m).*conj(array2D(s.phi_RIS_x(m),s.phi_RIS_y(m),cfg.RIS_size_x,cfg.RIS_size_y));
    for k=1:cfg.User_k
        response=s.alpha(k,m)*(s.Omega(:,:,m).'*array2D(s.th_RIS_x(k,m),s.th_RIS_y(k,m),cfg.RIS_size_x,cfg.RIS_size_y));
        s.rho(:,m)=s.rho(:,m)+s.X(:,k).*response;
        s.H(m,k,:)=response(cfg.T_p+1:end);
    end
end
s.Y_los=s.A_BS*s.rho.'; s.Y_noisefree=s.Y_los;
losPower=mean(abs(s.Y_los(:)).^2);
if strcmp(cfg.noise_mode,'fixed'), s.noise_power=0.001*10^(cfg.noise_dbm/10);
else, s.noise_power=losPower/10^(cfg.snr_db/10); end
s.actual_snr_db=10*log10(losPower/s.noise_power);
r=RandStream('mt19937ar','Seed',seed+6000000);
s.noise=sqrt(s.noise_power/2)*(randn(r,size(s.Y_los))+1j*randn(r,size(s.Y_los)));
[s.Y_nlos,s.nlos]=synthesize_nlos(cfg,s,seed);
s.interference_power=mean(abs(s.Y_nlos(:)).^2);
s.effective_noise_power=s.noise_power+s.interference_power;
s.Y=s.Y_los+s.Y_nlos+s.noise;
% Only the known TOTAL interference power enters the estimator covariance.
% Its spatial structure is deliberately unmodelled, as in the accepted setup.
end
