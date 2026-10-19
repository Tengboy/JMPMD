function ref=reference_g(cfg,scene)
User_k=cfg.User_k; RIS_m=cfg.RIS_m; T=cfg.T; T_d=cfg.T-cfg.T_p; pos_var=cfg.pos_var;
noise_power=scene.noise_power; trial=1; idx=1; pos_bound=0; th_bound=0; X_bound=NaN;
User_pos=scene.User_pos;
alpha=scene.alpha;
th_RIS_x=scene.th_RIS_x;
th_RIS_y=scene.th_RIS_y;
A_BS=scene.A_BS;
RIS_pos=scene.RIS_pos;
RIS_dir_x=scene.RIS_dir_x;
RIS_dir_y=scene.RIS_dir_y;
RIS_para=scene.RIS_para;
X=scene.X;
Omega=scene.Omega;
        %% Cramer Rao bound calculation
        % Position-prior precision: every user has covariance pos_var*I.
        G22 = [1/pos_var * eye(3),zeros(3,2*RIS_m);zeros(2*RIS_m,3),zeros(2*RIS_m,2*RIS_m)];
        for ind_u = 2:User_k
            G22 = blkdiag(G22,[1/pos_var * eye(3),zeros(3,2*RIS_m);zeros(2*RIS_m,3),zeros(2*RIS_m,2*RIS_m)]);
        end
        FIM_para_pos_pri = G22;
        FIM_para_pos_pri = blkdiag(FIM_para_pos_pri,2*eye(2*User_k*T_d));
        input.FIM_para_pos_pri = FIM_para_pos_pri;
        input.UE_location_r = User_pos;
        input.rho = alpha;
        input.noise_C_inv = (A_BS'*A_BS)/noise_power;
        input.th_RIS_x = th_RIS_x;
        input.th_RIS_y = th_RIS_y;
        input.X = X;
        Trans_mat1 = zeros((3+2*RIS_m)*User_k,4*RIS_m*User_k);
        %         Trans_mat1 = zeros((3+1*RIS_m)*User_k,3*RIS_m*User_k);
        for ind_r =1:RIS_m
            Trans_mat_temp = [];
            for ind_u = 1:User_k
                e_RU = (squeeze(User_pos(:,ind_u))-RIS_pos(:,ind_r))/norm(squeeze(User_pos(:,ind_u))-RIS_pos(:,ind_r));
                t1 = pi * (RIS_dir_x(:,ind_r) - e_RU.'*RIS_dir_x(:,ind_r)*e_RU)/norm(squeeze(User_pos(:,ind_u))-RIS_pos(:,ind_r));
                t2 = pi * (RIS_dir_y(:,ind_r) - e_RU.'*RIS_dir_y(:,ind_r)*e_RU)/norm(squeeze(User_pos(:,ind_u))-RIS_pos(:,ind_r));
                Trans_mat_temp2 = [t1,t2];
                Trans_mat_temp3 = [zeros((ind_r-1)*2,2);eye(2);zeros((RIS_m-ind_r)*2,2)];
                %                 Trans_mat_temp3 = [zeros((ind_r-1),1);eye(1);zeros((RIS_m-ind_r),1)];
                Trans_mat_temp = blkdiag(Trans_mat_temp,Trans_mat_temp2,Trans_mat_temp3);
            end
            Trans_mat1(:,(ind_r-1)*4*User_k+1:ind_r*4*User_k) = Trans_mat_temp;
            %             Trans_mat1(:,(ind_r-1)*3*User_k+1:ind_r*3*User_k) = Trans_mat_temp;
        end
        Trans_mat2 = eye(2*User_k*T_d);
        Trans_mat = blkdiag(Trans_mat1,Trans_mat2);
        
        input.Trans_mat = Trans_mat;
        input.Omega1 = Omega;
        input.RIS_para = RIS_para;
        input.T = T;
        input.T_d = T_d;
        [pos_bound(trial,idx),X_bound(trial,idx), FIM_para_pos,CRB]= Bound_cal_n_X3(input);
        % Only angle/gain columns are needed; avoid a full data-sized projection.
        parameterMap = sparse(Trans_mat(:,1:4*RIS_m*User_k));
        CRB_para_diag = sum(conj(parameterMap).*(CRB*parameterMap),1).';
        for ind_r = 1:RIS_m
            for ind_u = 1:User_k
                th_bound(trial,idx) = th_bound(trial,idx) +...
                    (CRB_para_diag((ind_r-1)*(4*User_k)+(ind_u-1)*4+1) + CRB_para_diag((ind_r-1)*(4*User_k)+(ind_u-1)*4+2));
            end
        end
        th_bound(trial,idx) = th_bound(trial,idx)/2/(RIS_m*User_k);
        
scale=1./sqrt(abs(diag(FIM_para_pos)));
scaledF=scale.*FIM_para_pos.*scale.';
scaledC=CRB./scale./scale.';
ref=struct('position_bound',pos_bound,'theta_bound',th_bound,'data_bound',X_bound, ...
    'fim_rcond',rcond(FIM_para_pos),'scaled_fim_rcond',rcond(scaledF), ...
    'inverse_residual',norm(scaledF*scaledC-eye(size(scaledF)),inf)/max(1,norm(scaledF,inf)*norm(scaledC,inf)));
ref.finite=all(isfinite(CRB(:)))&&all(real(diag(CRB))>0);
end
