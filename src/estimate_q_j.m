function result=estimate_q_j(cfg,scene,recordTrace)
% Adapted from the original JMPMD experiment Q_main_code_7_20_iter.m.
% Numerical update order is retained; see metadata/source_provenance.json.
if nargin<3, recordTrace=false; end
started=tic;
angle_solver=@Heuristic3;
% QPSK joint position / data experiment.
% theta variables use pi * direction cosine (spatial-frequency radians).
% Configuration/model choices are documented in docs/METHODS.md.
% exploiting the variational method to joint estimation line spectrum and transmit data

% Keep numerical warnings visible (especially singular/ill-conditioned FIMs).
% path setting
rng(scene.seed+10000000,'twister');
% Compare these modes with identical seeds and all other settings unchanged.
use_strict_moments = cfg.strict;
stop_tolerance = cfg.stop_tolerance;
stop_stable_iterations = cfg.stop_stable_iterations;
if use_strict_moments, result_suffix='strict'; else, result_suffix='approx'; end

%% phisical implement
f_c = cfg.f_c;
lambda = 3e8/f_c;

Tx_power_dB = cfg.Tx_power_dB;

SNR = cfg.snr_db;
User_k = cfg.User_k;
RIS_m = cfg.RIS_m;
mont = 1;
pos_var = cfg.pos_var;
counter_max = 1;
counter_max0 = cfg.max_iterations;
T = cfg.T;
T_p = cfg.T_p;
T_d = T - T_p;       % data timeslot num
constellation = [1+1j,-1+1j,-1-1j,1-1j]/sqrt(2);    % QPSK

X_state = constellation(Constell_Mat(User_k,length(constellation)));

%% BS parameter
BS_pos = [-20,0,0].';   % BS location
BS_dir = [0,1,0].';     % antenna direction of BS
BS_size = cfg.BS_size;

%% RIS parameter
RIS_size_x = cfg.RIS_size_x;
RIS_size_y = cfg.RIS_size_y;
RIS_size = RIS_size_x * RIS_size_y;
damp_a = 0.1;

RIS_pos = [0,20,10;
    0,-20,10;
    20,0,5].';
RIS_pos = RIS_pos(:,1:RIS_m);

RIS_dir_x = [1,0,0;
    -1,0,0;
    0,1,0].';
RIS_dir_x = RIS_dir_x(:,1:RIS_m);

RIS_dir_y = [0,0,1;
    0,0,1;
    0,0,1].';
RIS_dir_y = RIS_dir_y(:,1:RIS_m);

RIS_dir_z = [0,-1,0;
    0,1,0;
    -1,0,0].';
RIS_dir_z = RIS_dir_z(:,1:RIS_m);

% anchor parameters used in taylor localization
RIS_para.RIS_dir_x = RIS_dir_x;
RIS_para.RIS_dir_y = RIS_dir_y;
RIS_para.RIS_size_x = RIS_size_x;
RIS_para.RIS_size_y = RIS_size_y;
RIS_para.RIS_pos = RIS_pos;

%% UE pos
User_pos_c = [-5,0,3.5;
    10,10,1;
    10,-10,1].';            % center location of users

%% preprocessing
phi_RIS_x = pi * diag(RIS_dir_x.' * (BS_pos - RIS_pos))./vecnorm(BS_pos - RIS_pos).';   % AoD of RIS
phi_RIS_y = pi * diag(RIS_dir_y.' * (BS_pos - RIS_pos))./vecnorm(BS_pos - RIS_pos).';   % AoD of RIS
th_BS = pi * (BS_dir.' * (RIS_pos - BS_pos)).' ./vecnorm(BS_pos - RIS_pos).';           % AoA of BS
d_B = vecnorm(BS_pos - RIS_pos).';                                                      % distance of BS and RIS
A_BS = array(th_BS,BS_size);

pos_MSE = zeros(mont,length(SNR));
pos_MSE_iter_sum = zeros(counter_max0,length(SNR));
th_MSE = zeros(mont,length(SNR));
th_MSE_iter_sum = zeros(counter_max0,length(SNR));
G_MSE_iter_sum = zeros(counter_max0,length(SNR));
th_f_MSE = zeros(mont,length(SNR));
% Historical data_BER* names store SYMBOL error rates (SER), not bit error rates.
data_BER1 = zeros(mont,length(SNR));
data_BER2 = zeros(mont,length(SNR));
data_BER3 = zeros(mont,length(SNR)); % 完美

LM_MSE = zeros(counter_max0,1);
LM_MSE_e = zeros(counter_max0,1);

VR_MSE = zeros(counter_max0,1);
VR_MSE_end = zeros(mont,length(SNR));
VR_MSE_e = zeros(counter_max0,1);

VR_PRI_MSE = zeros(counter_max0,1);
VR_PRI_MSE_e = zeros(counter_max0,1);

LM_PRI_MSE = zeros(counter_max0,1);
LM_PRI_MSE_e = zeros(counter_max0,1);

pos_bound = zeros(mont,length(SNR));
th_bound = zeros(mont,length(SNR));

% Post-update metrics: G MSE, G NMSE, data-symbol soft MSE,
% spatial-frequency MSE, position MSE. Pilot-only traces use pilot G only.
convergence_metric_names = {'G_MSE','G_NMSE','X_MSE','theta_MSE','position_MSE'};
convergence_sum = zeros(counter_max0,5,length(SNR));
convergence_count = zeros(counter_max0,5,length(SNR));
convergence_active_count = zeros(counter_max0,length(SNR));
iteration_count = zeros(mont,length(SNR));
iteration_converged = false(mont,length(SNR));
localization_calls = zeros(mont,length(SNR));
localization_failures = zeros(mont,length(SNR));
localization_failure_reasons = zeros(5,length(SNR));

for idx = 1:length(SNR)
    Tx_power = 0.001 * 10^(Tx_power_dB/10);     % transmit power (W)
    discoverage = [];
    for trial = 1:mont

        %% Shared observation: only pilot/data availability differs by method.
        iter_flag0=1;
        User_pos=scene.User_pos;
        th_RIS_x=scene.th_RIS_x;
        th_RIS_y=scene.th_RIS_y;
        th_RIS=scene.th_RIS;
        d_U=scene.d_U;
        alpha=scene.alpha;
        X=scene.X;
        X_p=scene.X_p;
        X_d=scene.X_d;
        Omega=scene.Omega;
        rho=scene.rho;
        H=scene.H;
        Y=scene.Y;
        Y_noisefree=scene.Y_noisefree;
        data_ind_p=scene.data_ind_p;
        data_ind_d=scene.data_ind_d;
        noise_power=scene.effective_noise_power;
        R_W=zeros(RIS_size,RIS_size,RIS_m);
        for R_idx=1:RIS_m
            R_W(:,:,R_idx)=conj(Omega(:,:,R_idx))*Omega(:,:,R_idx).';
        end
        sR_W=ris_lag_sum(R_W,RIS_size_x,RIS_size_y);
        %% estimate rho (LS estimation)
        
        G_B_pos_C = zeros(T,RIS_m);
        G_B_ext_m = zeros(T,RIS_m);
        G_B_ext_C = zeros(T,RIS_m);
        G_A_pos_m = zeros(T,RIS_m);
        G_LM_C_r = zeros(T,RIS_m,RIS_m);
        G_A_pos_C = zeros(T,RIS_m);
        
        G_A_ext_m = ((A_BS'*A_BS)\A_BS' * Y).';
        G_A_ext_C = diag(diag(noise_power*inv(A_BS'*A_BS)));
        G_B_ext_m_ini = G_A_ext_m;
        G_B_ext_C_ini = ones(T,RIS_m)*G_A_ext_C;

        %% initial hyperparameters nu and tau
        nu = diag(G_A_ext_C);
        tau = 10000*ones(RIS_m,1);
        
        %% initial X_d as 0 (random sampling ?) (flat prior) ??
        q_x_p = zeros(4,T_p,User_k);
        for T_idx = 1:T_p
            for U_idx = 1:User_k
                q_x_p(data_ind_p(T_idx,U_idx),T_idx,U_idx) = 1;
            end
        end
        q_x_d = zeros(4,T_d,User_k);
        data_ind_d_ini = randi(4,[T_d,User_k]);
        for T_idx = 1:T_d
            for U_idx = 1:User_k
                q_x_d(:,T_idx,U_idx)=1/length(constellation);
            end
        end
        q_x = cat(2,q_x_p,q_x_d);
        
        X_est_ini = zeros(T,User_k);
        X_est_var_ini = zeros(T,User_k,User_k);
        for U_idx = 1:User_k
            X_est_ini(:,U_idx) = (constellation * q_x(:,:,U_idx)).';
            X_est_var_ini(:,U_idx,U_idx) = sum(abs(X_est_ini(:,U_idx) - constellation).^2 .* q_x(:,:,U_idx).',2);
        end

        % initial theta and alpha using pilot
        % initial theta method 1, using the prior position information
        th_est = zeros(User_k,RIS_m,2);
        th_est_ka = zeros(User_k,RIS_m,2);
        th_pos = zeros(User_k,RIS_m,2);
        th_ka_pos = zeros(User_k,RIS_m,2);
        A_th_x = zeros(RIS_size_x,User_k,RIS_m);
        A_th_y = zeros(RIS_size_y,User_k,RIS_m);
        A_th_x_exd = zeros(2*RIS_size_x-1,User_k,RIS_m);
        A_th_y_exd = zeros(2*RIS_size_y-1,User_k,RIS_m);
        A_th_exd = zeros((2*RIS_size_y-1)*(2*RIS_size_x-1),User_k,RIS_m);
        A_th = zeros(RIS_size,User_k,RIS_m);
        A_a = zeros(T,User_k,RIS_m);
        for R_idx = 1:RIS_m
            for U_idx = 1:User_k
                UE_cov_pri = pos_var * eye(3);
                d_U_pri = norm(User_pos_c(:,U_idx) - RIS_pos(:,R_idx));
                e3 = (User_pos_c(:,U_idx) - RIS_pos(:,R_idx))/d_U_pri;
                e1 = cross(e3,cross(RIS_dir_x(:,R_idx),e3));
                e1 = e1/norm(e1);
                e2 = cross(e3,cross(RIS_dir_y(:,R_idx),e3));
                e2 = e2/norm(e2);
                th_est(U_idx,R_idx,1) = pi * e3.'*RIS_dir_x(:,R_idx);
                th_est(U_idx,R_idx,2) = pi * e3.'*RIS_dir_y(:,R_idx);
                th_est_ka(U_idx,R_idx,1) = d_U_pri^2/(e1.' * UE_cov_pri * e1)/(pi^2 - th_est(U_idx,R_idx,1)^2);
                th_est_ka(U_idx,R_idx,2) = d_U_pri^2/(e2.' * UE_cov_pri * e2)/(pi^2 - th_est(U_idx,R_idx,2)^2);
                %                 th_est(U_idx,R_idx,1) = th_est(U_idx,R_idx,1) - phi_RIS_x(R_idx);
                %                 th_est(U_idx,R_idx,2) = th_est(U_idx,R_idx,2) - phi_RIS_y(R_idx);
                if isfield(cfg,'initial_theta')
                    th_est(U_idx,R_idx,:)=cfg.initial_theta(U_idx,R_idx,:);
                end
                A_th_x(:,U_idx,R_idx) = array_B(th_est(U_idx,R_idx,1),th_est_ka(U_idx,R_idx,1),RIS_size_x);
                A_th_y(:,U_idx,R_idx) = array_B(th_est(U_idx,R_idx,2),th_est_ka(U_idx,R_idx,2),RIS_size_y);
                A_th(:,U_idx,R_idx) = array2D_B(th_est(U_idx,R_idx,1),th_est(U_idx,R_idx,2),...
                    th_est_ka(U_idx,R_idx,1),th_est_ka(U_idx,R_idx,2),RIS_size_x,RIS_size_y);
                A_a(:,U_idx,R_idx) = X_est_ini(:,U_idx) .* (Omega(:,:,R_idx).'*A_th(:,U_idx,R_idx));
            end
            A_th_x_exd(:,:,R_idx) = [flipud(A_th_x(2:end,:,R_idx))'.';A_th_x(:,:,R_idx)];
            A_th_y_exd(:,:,R_idx) = [flipud(A_th_y(2:end,:,R_idx))'.';A_th_y(:,:,R_idx)];
            A_th_exd(:,:,R_idx) = KR(A_th_x_exd(:,:,R_idx),A_th_y_exd(:,:,R_idx));
        end
        
        % initial alpha with LMMSE
        alpha_est = zeros(User_k,RIS_m);
        alpha_est_cov = zeros(User_k,User_k,RIS_m);
        for R_idx = 1:RIS_m
            J_temp = A_a(:,:,R_idx)' * A_a(:,:,R_idx);
            if use_strict_moments
                % The first strict gain distribution uses known pilots only.
                % Unknown symbols have zero initial mean; h already uses pilots.
                first_slots=1:T;
                if cfg.gain_pilot_init,first_slots=1:T_p;end
                [psiMean,psiPower] = vb_steering_moments(Omega(:,first_slots,R_idx), ...
                    A_th_x(:,:,R_idx),A_th_y(:,:,R_idx));
                J_temp = vb_gain_gram(X_est_ini(first_slots,:),permute(X_est_var_ini(first_slots,:,:),[2 3 1]),psiMean,psiPower);
            end
            h_temp = A_a(:,:,R_idx)' * G_A_ext_m(:,R_idx);
            alpha_est_cov(:,:,R_idx) = inv(J_temp/nu(R_idx)+eye(User_k)/tau(R_idx));
            alpha_est(:,R_idx) = alpha_est_cov(:,:,R_idx) * h_temp/nu(R_idx);
        end
        
        %% initial X_d given q(theta) and q(alpha)

        H_e = zeros(RIS_m,User_k,T_d);
        H_e_b = zeros(RIS_m,User_k,T_d);
        X_est = X_est_ini;
        X_est_var = zeros(T,User_k,User_k);

        p_like = zeros(length(X_state),T_d);
        if use_strict_moments
            HH_initial = vb_data_gram(Omega(:,T_p+1:end,:),A_th_x,A_th_y, ...
                alpha_est,alpha_est_cov,inv(G_A_ext_C));
        end
        for T_idx = 1:T_d
            for R_idx = 1:RIS_m
                for U_idx = 1:User_k
                    H_e(R_idx,U_idx,T_idx) =alpha_est(U_idx,R_idx) * Omega(:,T_p+T_idx,R_idx).' * A_th(:,U_idx,R_idx);
                end
            end
            %             H_nmse(T_idx) = norm(H_e(:,:,T_idx)-H(:,:,T_idx),'fro')/norm(H(:,:,T_idx),'fro')
            
            H_e_temp = H_e(:,:,T_idx);
            rho_est_temp = G_A_ext_m(T_p+T_idx,:).';

            residual = rho_est_temp - H_e_temp*X_state.';
            log_p_like = -abs(sum(conj(residual).*(G_A_ext_C\residual),1)).';
            if use_strict_moments
                linearTerm = H_e_temp'*(G_A_ext_C\rho_est_temp);
                log_p_like = -real(sum((conj(X_state)*HH_initial(:,:,T_idx)).*X_state,2)) ...
                    + 2*real(conj(X_state)*linearTerm);
            end
            log_p_like = log_p_like - max(log_p_like); % 防止exp的数值问题
            p_like(:,T_idx) = exp(log_p_like);
            p_like(:,T_idx) = p_like(:,T_idx)/sum(p_like(:,T_idx));    % 归一化
            [X_est(T_p+T_idx,:),X_est_var(T_p+T_idx,:,:)] = ...
                discrete_posterior_moments(X_state,p_like(:,T_idx));
        end
        
        for R_idx = 1:RIS_m
            for U_idx = 1:User_k
                A_a(:,U_idx,R_idx) = X_est(:,U_idx) .* (Omega(:,:,R_idx).'*A_th(:,U_idx,R_idx));
            end
        end
        % initial alpha with LMMSE
        for R_idx = 1:RIS_m
            J_temp = A_a(:,:,R_idx)' * A_a(:,:,R_idx);
            if use_strict_moments
                [psiMean,psiPower] = vb_steering_moments(Omega(:,1:T,R_idx), ...
                    A_th_x(:,:,R_idx),A_th_y(:,:,R_idx));
                J_temp = vb_gain_gram(X_est(1:T,:),permute(X_est_var(1:T,:,:),[2 3 1]),psiMean,psiPower);
            end
            h_temp = A_a(:,:,R_idx)' * G_A_ext_m(:,R_idx);
            alpha_est_cov(:,:,R_idx) = inv(J_temp/nu(R_idx)+eye(User_k)/tau(R_idx));
            alpha_est(:,R_idx) = alpha_est_cov(:,:,R_idx) * h_temp/nu(R_idx);
        end
        %% message iteration
        
        mse_out_o = 0;
        th_pri = th_est;
        th_ka_pri = th_est_ka;
        th_pri_eta = th_ka_pri.*exp(1j*th_pri);
        
        counter0 = 0;
        mse_out = zeros(counter_max0,1);
        th_MSE_iter = zeros(counter_max0,1);
        pos_MSE_iter = zeros(counter_max0,1);
        X_MSE = zeros(counter_max0,1);
%         VR_PRI_MSE(counter0) = sum(abs(G_A_ext_m - rho).^2,'all');
%         VR_PRI_MSE_e(counter0) = sum(G_A_ext_C,'all')*50;
        stable_iterations = 0;
        trace_metrics = nan(counter_max0,5);
        previous_state.G = zeros(T,RIS_m);
        for R_idx=1:RIS_m
            previous_state.G(:,R_idx) = A_a(:,:,R_idx)*alpha_est(:,R_idx);
        end
        previous_state.X = X_est(1:T,:);
        previous_state.position = User_pos_c.';
        completed_iterations = 0;
        while iter_flag0
            counter0 = counter0+1;
            iter_flag = 1;
            counter = 1;
            mse = zeros(counter_max,1);
            alpha_MSE = zeros(counter_max,1);
            
            X_MSE(1) = sum(abs(X_est(T_p+1:end,:)-X(T_p+1:end,:)).^2,'all')/(T_d*User_k);
            
            rho_denoise = zeros(T,RIS_m);
            for R_idx = 1:RIS_m
                rho_denoise(:,R_idx) = A_a(:,:,R_idx)*alpha_est(:,R_idx);
            end
            mse(counter) = norm(rho(1:T_p,:) - rho_denoise(1:T_p,:),'fro')^2/norm(rho(1:T_p,:),'fro')^2;
            %% update q(alpha) for all RIS and all user and EM for nu, tau
            for R_idx = 1:RIS_m
                J_temp = A_a(:,:,R_idx)' * A_a(:,:,R_idx);
                sR_W_temp = sR_W(:,:,R_idx);
                J_temp = J_temp - diag(diag(J_temp)) + diag((real(sR_W_temp(:).'*A_th_exd(:,:,R_idx))).');
                if use_strict_moments
                    [psiMean,psiPower] = vb_steering_moments(Omega(:,1:T,R_idx), ...
                        A_th_x(:,:,R_idx),A_th_y(:,:,R_idx));
                    J_temp = vb_gain_gram(X_est(1:T,:),permute(X_est_var(1:T,:,:),[2 3 1]),psiMean,psiPower);
                end
                h_temp = A_a(:,:,R_idx)' * G_A_ext_m(:,R_idx);
                alpha_est_cov(:,:,R_idx) = inv(J_temp/nu(R_idx)+eye(User_k)/tau(R_idx));
                alpha_est(:,R_idx) = alpha_est_cov(:,:,R_idx) * h_temp/nu(R_idx);
                y2_p = G_A_ext_m(:,R_idx)'* G_A_ext_m(:,R_idx);
                %                     nu(R_idx)  = real(y2_p - 2*real(h_temp'*alpha_est(:,R_idx)) + alpha_est(:,R_idx)'*J_temp*alpha_est(:,R_idx) + trace(J_temp*alpha_est_cov(:,:,R_idx)))/T;
                tau(R_idx) = real(alpha_est(:,R_idx)'*alpha_est(:,R_idx)+trace(alpha_est_cov(:,:,R_idx)) )/User_k;
            end
            alpha_MSE(counter) = norm(alpha_est - alpha,'fro');
            %% update q(theta) for all RIS and all user
            for R_idx = 1:RIS_m
                for U_idx = 1:User_k
                    square_exp = real((alpha_est_cov(U_idx,U_idx,R_idx)+abs(alpha_est(U_idx,R_idx)).^2));
                    square_exp_i = conj(alpha_est_cov([1:U_idx-1 U_idx+1:end],U_idx,R_idx) + alpha_est([1:U_idx-1 U_idx+1:end],R_idx)*alpha_est(U_idx,R_idx)');
                    eta_th_cov_mat = zeros(User_k,RIS_size);
                    for U_idx2 = 1:User_k
                        projectedSteering = A_th(:,U_idx2,R_idx)' * conj(Omega(:,:,R_idx));
                        eta_th_cov_mat(U_idx2,:) = (projectedSteering .* ...
                            reshape(X_est_var(:,U_idx,U_idx2),1,[])) * Omega(:,:,R_idx).';
                    end
                    matchedResidual = G_A_ext_m(:,R_idx)'*alpha_est(U_idx,R_idx) ...
                        - square_exp_i.'*A_a(:,[1:U_idx-1,U_idx+1:User_k],R_idx)';
                    eta_th = 2/nu(R_idx)*((matchedResidual .* X_est(:,U_idx).')*Omega(:,:,R_idx).' ...
                        - square_exp_i.'*eta_th_cov_mat([1:U_idx-1 U_idx+1:end],:))';
                    Mat_eta_th = reshape(eta_th'.',RIS_size_y,RIS_size_x);
                    eta_th_x = (A_th_y(:,U_idx,R_idx).' * Mat_eta_th)';
                    eta_th_y = (Mat_eta_th*A_th_x(:,U_idx,R_idx))'.';
                    % part 2
                    sR_temp_x = sR_W(:,:,R_idx).'*A_th_y_exd(:,U_idx,R_idx);
                    sR_temp_x = sR_temp_x(RIS_size_x:end);
                    sR_temp_x(2:end) = 2 * sR_temp_x(2:end);
                    sR_temp_x = sR_temp_x'.';
                    sR_temp_y = sR_W(:,:,R_idx)*A_th_x_exd(:,U_idx,R_idx);
                    sR_temp_y = sR_temp_y(RIS_size_y:end);
                    sR_temp_y(2:end) = 2 * sR_temp_y(2:end);
                    sR_temp_y = sR_temp_y'.';
                    eta_th_x = eta_th_x - 1/nu(R_idx) * square_exp * sR_temp_x;
                    eta_th_y = eta_th_y - 1/nu(R_idx) * square_exp * sR_temp_y;
                    ind_x = find(abs(eta_th_x)>0)-1;
                    ind_y = find(abs(eta_th_y)>0)-1;
                    %                         [A_th_x(:,U_idx,R_idx), th_est(U_idx,R_idx,1), th_est_ka(U_idx,R_idx,1)] = angle_solver(eta_th_x,ind_x,th_pri(U_idx,R_idx,1),th_ka_pri(U_idx,R_idx,1),th_est(U_idx,R_idx,1));
                    %                         [A_th_y(:,U_idx,R_idx), th_est(U_idx,R_idx,2), th_est_ka(U_idx,R_idx,2)] = angle_solver(eta_th_y,ind_y,th_pri(U_idx,R_idx,2),th_ka_pri(U_idx,R_idx,2),th_est(U_idx,R_idx,2));
                    [A_th_x(:,U_idx,R_idx), th_est(U_idx,R_idx,1), th_est_ka(U_idx,R_idx,1)] = angle_solver(eta_th_x,ind_x,0,0,th_est(U_idx,R_idx,1));
                    [A_th_y(:,U_idx,R_idx), th_est(U_idx,R_idx,2), th_est_ka(U_idx,R_idx,2)] = angle_solver(eta_th_y,ind_y,0,0,th_est(U_idx,R_idx,2));
%                     A_th(:,U_idx,R_idx) = KR(A_th_x(:,U_idx,R_idx),A_th_y(:,U_idx,R_idx));
%                     A_a(:,U_idx,R_idx) = X_est(:,U_idx) .* (Omega(:,:,R_idx).'*A_th(:,U_idx,R_idx));
                end
%                 A_th_x_exd(:,:,R_idx) = [flipud(A_th_x(2:end,:,R_idx))'.';A_th_x(:,:,R_idx)];
%                 A_th_y_exd(:,:,R_idx) = [flipud(A_th_y(2:end,:,R_idx))'.';A_th_y(:,:,R_idx)];
%                 A_th_exd(:,:,R_idx) = KR(A_th_x_exd(:,:,R_idx),A_th_y_exd(:,:,R_idx));
            end
            
            %%

            th_ext = th_est;
            th_ka_ext = th_est_ka;
            
            %% position estimation (merge messages from theta to position)
            User_pos_est = zeros(User_k,3);
            User_pos_cov = zeros(3,3,User_k);
            for U_idx = 1:User_k
                AoA_est = squeeze(th_ext(U_idx,:,:));
                [User_pos_est(U_idx,:),User_pos_cov(:,:,U_idx),~,localization_status] = taylor_localization(AoA_est,RIS_para,squeeze(1./th_ka_ext(U_idx,:,:)),User_pos_c(:,U_idx));
                localization_calls(trial,idx) = localization_calls(trial,idx)+1;
                if ~localization_status.converged
                    localization_failures(trial,idx) = localization_failures(trial,idx)+1;
                    reason = localization_status.reason_code;
                    localization_failure_reasons(reason,idx) = localization_failure_reasons(reason,idx)+1;
                end
                % merge with pri
                User_pos_cov_n = inv(inv(User_pos_cov(:,:,U_idx)) + inv(UE_cov_pri));
                User_pos_est_n = User_pos_cov_n * (User_pos_cov(:,:,U_idx)\User_pos_est(U_idx,:).' + UE_cov_pri \ User_pos_c(:,U_idx));
                User_pos_est(U_idx,:) = User_pos_est_n;
                User_pos_cov(:,:,U_idx) = User_pos_cov_n;
                
                % update theta_est
                for R_idx = 1:RIS_m
                    d_U_pri = norm(User_pos_est(U_idx,:).' - RIS_pos(:,R_idx));
                    e3 = (User_pos_est(U_idx,:).' - RIS_pos(:,R_idx))/d_U_pri;
                    e1 = cross(e3,cross(RIS_dir_x(:,R_idx),e3));
                    e1 = e1/norm(e1);
                    e2 = cross(e3,cross(RIS_dir_y(:,R_idx),e3));
                    e2 = e2/norm(e2);
                    th_pos(U_idx,R_idx,1) = pi * e3.'*RIS_dir_x(:,R_idx);
                    th_pos(U_idx,R_idx,2) = pi * e3.'*RIS_dir_y(:,R_idx);
                    th_ka_pos(U_idx,R_idx,1) = d_U_pri^2/(e1.' * User_pos_cov(:,:,U_idx) * e1)/(pi^2 - th_pos(U_idx,R_idx,1)^2);
                    th_ka_pos(U_idx,R_idx,2) = d_U_pri^2/(e2.' * User_pos_cov(:,:,U_idx) * e2)/(pi^2 - th_pos(U_idx,R_idx,2)^2);

                    A_th_x(:,U_idx,R_idx) = array_B(th_pos(U_idx,R_idx,1),th_ka_pos(U_idx,R_idx,1),RIS_size_x);
                    A_th_y(:,U_idx,R_idx) = array_B(th_pos(U_idx,R_idx,2),th_ka_pos(U_idx,R_idx,2),RIS_size_y);
                    A_th(:,U_idx,R_idx) = KR(A_th_x(:,U_idx,R_idx),A_th_y(:,U_idx,R_idx));
                    A_a(:,U_idx,R_idx) = X_est(:,U_idx) .* (Omega(:,:,R_idx).'*A_th(:,U_idx,R_idx));
                    A_th_x_exd(:,U_idx,R_idx) = [flipud(A_th_x(2:end,U_idx,R_idx))'.';A_th_x(:,U_idx,R_idx)];
                    A_th_y_exd(:,U_idx,R_idx) = [flipud(A_th_y(2:end,U_idx,R_idx))'.';A_th_y(:,U_idx,R_idx)];
                    A_th_exd(:,U_idx,R_idx) = KR(A_th_x_exd(:,U_idx,R_idx),A_th_y_exd(:,U_idx,R_idx));
                end
            end
            th_MSE_iter(counter0) = sum((th_pos(:,:,1:2) - th_RIS(:,:,1:2)).^2,'all');
            pos_MSE_iter(counter0) = sum((User_pos_est.' - User_pos).^2,'all')/User_k;

            HH_mat = zeros(User_k,User_k,T_d);
            Hy_mat = zeros(User_k,T_d);
            
            H2_mat = zeros(RIS_m,User_k,T_d);
            for R_idx = 1:RIS_m
                for U_idx = 1:User_k
                    H_e_b(R_idx,U_idx,:) =alpha_est(U_idx,R_idx) * Omega(:,T_p+1:end,R_idx).' * array2D_B(th_pos(U_idx,R_idx,1),th_pos(U_idx,R_idx,2),...
                        th_ka_pos(U_idx,R_idx,1),th_ka_pos(U_idx,R_idx,2),RIS_size_x,RIS_size_y);
                    
                    % E[|w^T a|^2]; compute its diagonal without a Td-by-Td temporary.
                    meanPower = steering_power(Omega(:,T_p+1:end,R_idx), ...
                    array_B(th_pos(U_idx,R_idx,1),th_ka_pos(U_idx,R_idx,1),RIS_size_x), array_B(th_pos(U_idx,R_idx,2),th_ka_pos(U_idx,R_idx,2),RIS_size_y)).';
                    H2_mat(R_idx,U_idx,:) = real((alpha_est_cov(U_idx,U_idx,R_idx)+abs(alpha_est(U_idx,R_idx)).^2)) * meanPower;
                end
            end
            if use_strict_moments
                HH_strict = vb_data_gram(Omega(:,T_p+1:end,:),A_th_x,A_th_y, ...
                    alpha_est,alpha_est_cov,diag(1./nu));
            end
            for T_idx = 1:T_d
                HH_mat(:,:,T_idx) = H_e_b(:,:,T_idx)'*diag(1./nu)*H_e_b(:,:,T_idx);
                HH_mat(:,:,T_idx) = HH_mat(:,:,T_idx) - diag(diag(H_e_b(:,:,T_idx)'*diag(1./nu)*H_e_b(:,:,T_idx))) + diag((1./nu).'*H2_mat(:,:,T_idx));
                if use_strict_moments
                    HH_mat(:,:,T_idx) = HH_strict(:,:,T_idx);
                end

                Hy_mat(:,T_idx) = G_A_ext_m(T_p+T_idx,:)'.'*diag(1./nu)*H_e_b(:,:,T_idx);

                % Evaluate all joint QPSK states together; retain exact enumeration.
                log_p_like = -real(sum((conj(X_state)*HH_mat(:,:,T_idx)).*X_state,2)) ...
                    + 2*real(X_state*Hy_mat(:,T_idx));
                log_p_like = log_p_like - max(log_p_like); % 防止exp的数值问题
                p_like(:,T_idx) = exp(log_p_like);
                p_like(:,T_idx) = p_like(:,T_idx)/sum(p_like(:,T_idx));    % 归一化
                [X_est(T_p+T_idx,:),X_est_var(T_p+T_idx,:,:)] = ...
                    discrete_posterior_moments(X_state,p_like(:,T_idx));
            end
            
            for R_idx= 1:RIS_m
                for U_idx = 1:User_k
                    A_a(:,U_idx,R_idx) = X_est(:,U_idx) .* (Omega(:,:,R_idx).'*A_th(:,U_idx,R_idx));
                end
            end
            
            X_hat_hard = zeros(T_d,User_k);
            for T_idx = 1:T_d
                [~,max_idx]=max(p_like(:,T_idx));
                X_hat_hard(T_idx,:) = X_state(max_idx,:);
            end
            X_MSE(counter0) = sum(X_hat_hard~=X_d,'all')/(T_d*User_k);
            
            %%
            
            % Record AFTER gain, angle/position and data updates in this sweep.
            completed_iterations = completed_iterations+1;
            current_state.G = zeros(T,RIS_m);
            for R_idx=1:RIS_m
                current_state.G(:,R_idx) = A_a(:,:,R_idx)*alpha_est(:,R_idx);
            end
            current_state.X = X_est(1:T,:);
            current_state.position = User_pos_est;
            [stable_iterations,relative_change] = iteration_status(current_state,previous_state, ...
                stable_iterations,stop_tolerance);
            previous_state = current_state;
            g_error = norm(current_state.G-rho(1:T,:),'fro')^2;
            trace_metrics(completed_iterations,1:2) = [g_error/(T*RIS_m), ...
                g_error/max(norm(rho(1:T,:),'fro')^2,realmin)];
            trace_metrics(completed_iterations,3) = mean(abs(X_est(T_p+1:end,:)-X(T_p+1:end,:)).^2,'all');
            trace_metrics(completed_iterations,4) = mean((th_pos(:,:,1:2)-th_RIS(:,:,1:2)).^2,'all');
            trace_metrics(completed_iterations,5) = norm(User_pos_est.'-User_pos,'fro')^2/User_k;
            convergence_active_count(completed_iterations,idx) = convergence_active_count(completed_iterations,idx)+1;
            if stable_iterations>=stop_stable_iterations || completed_iterations>=counter_max0
                iter_flag0 = 0;
                iteration_count(trial,idx) = completed_iterations;
                iteration_converged(trial,idx) = stable_iterations>=stop_stable_iterations;
                % Hold the final value for equal trial counts in convergence curves.
                trace_metrics(completed_iterations+1:end,:) = repmat(trace_metrics(completed_iterations,:), ...
                    counter_max0-completed_iterations,1);
                % Only structurally unavailable metrics are omitted. Failed
                % numerical trials propagate NaN/Inf instead of disappearing.
                metric_available = true(1,5);
                metric_available(3) = true;
                metric_available(5) = true;
                valid_metric = repmat(metric_available,counter_max0,1);
                metrics_to_sum = trace_metrics; metrics_to_sum(~valid_metric)=0;
                convergence_sum(:,:,idx) = convergence_sum(:,:,idx)+metrics_to_sum;
                convergence_count(:,:,idx) = convergence_count(:,:,idx)+valid_metric;
            end
            mse_out = trace_metrics(:,2);
            X_MSE = trace_metrics(:,3);
            th_MSE_iter = trace_metrics(:,4);
            pos_MSE_iter = trace_metrics(:,5);
        end
        %% hard decision
        X_hat_hard = zeros(T_d,User_k);
        for T_idx = 1:T_d
            [~,max_idx]=max(p_like(:,T_idx));
            X_hat_hard(T_idx,:) = X_state(max_idx,:);
        end
        data_BER1(trial,idx) = sum(X_hat_hard~=X_d,'all')/(T_d*User_k);
        %% record results
        th_MSE(trial,idx) = sum((th_pos(:,:,1:2) - th_RIS(:,:,1:2)).^2,'all')/2/(User_k*RIS_m);
        pos_MSE(trial,idx) = sum((User_pos_est.' - User_pos).^2,'all')/User_k;
        pos_MSE_iter_sum(:,idx) = pos_MSE_iter_sum(:,idx) + pos_MSE_iter;
        th_MSE_iter_sum(:,idx) = th_MSE_iter_sum(:,idx) + th_MSE_iter;
        G_MSE_iter_sum(:,idx) = G_MSE_iter_sum(:,idx) + mse_out;
    end
end

data_BER_ave1 = mean(data_BER1,1);
data_BER_ave2 = mean(data_BER2,1);
data_BER_ave3 = mean(data_BER3,1);
th_MSE_ave = mean(th_MSE,1);
th_f_MSE_ave = mean(th_f_MSE,1);
pos_MSE_ave = mean(pos_MSE,1);
th_bound_ave = mean(th_bound,1);
pos_bound_ave = mean(pos_bound,1);
convergence_mean = convergence_sum./convergence_count;
localization_call_failure_rate = sum(localization_failures,1)./sum(localization_calls,1);
localization_trial_failure_rate = mean(localization_failures>0,1);
iteration_convergence_rate = mean(iteration_converged,1);
% Report both the number of calls and trials; do not silently remove failures.


result.theta_mse=th_MSE(1);
result.theta_post_mse=th_MSE(1);
result.position_user_error=sum((User_pos_est.'-User_pos).^2,1);
result.position_mse=mean(result.position_user_error);
result.data_error=data_BER1(1);
result.iterations=iteration_count(1); result.converged=iteration_converged(1);
result.cap=result.iterations>=cfg.max_iterations;
result.last_change=relative_change;
result.localization_calls=localization_calls(1);
result.localization_failures=localization_failures(1);
result.localization_reasons=localization_failure_reasons(:,1).';
result.finite=all(isfinite([result.theta_mse result.position_user_error result.data_error]));
result.trace=[];
if recordTrace, result.trace=trace_metrics; end
result.elapsed=toc(started);
end
