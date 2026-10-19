% Joint position / complex-symbol information bound.
function [obj1,obj2,FIM_para_pos,CRB] = Bound_cal_n_X3(input)
User_N = size(input.UE_location_r,2);
RIS_K = size(input.RIS_para.RIS_pos,2);
Nx = input.RIS_para.RIS_size_x;
Ny = input.RIS_para.RIS_size_y;
T = input.T;
T_d = input.T_d;
T_p = T-T_d;
eta_size = 4*User_N*RIS_K + 2*User_N*T_d;

% Channel coordinates are RIS-major [theta_x,theta_y,|rho|,arg(rho)].
% Data coordinates are user-major [real(x_d);imag(x_d)]. Each data
% derivative occupies only its own time row; retain that sparsity.
partial_mat = cell(RIS_K,1);
for ind_r = 1:RIS_K
    partial_mat{ind_r} = spalloc(T,eta_size,4*User_N*T+2*User_N*T_d);
    Wt = input.Omega1(:,:,ind_r).';
    for ind_u = 1:User_N
        rho = input.rho(ind_u,ind_r);
        ax = linespec(input.th_RIS_x(ind_u,ind_r),Nx);
        ay = linespec(input.th_RIS_y(ind_u,ind_r),Ny);
        a = KR(ax,ay);
        dax = KR((1j*(0:Nx-1).').*ax,ay);
        day = KR(ax,(1j*(0:Ny-1).').*ay);
        projected = Wt*[dax,day,a];
        weighted = input.X(:,ind_u).*projected;
        first = (ind_r-1)*4*User_N + (ind_u-1)*4;
        partial_mat{ind_r}(:,first+(1:4)) = ...
            [rho*weighted(:,1:2), (rho/abs(rho))*weighted(:,3), 1j*rho*weighted(:,3)];

        first_data = 4*User_N*RIS_K + (ind_u-1)*2*T_d;
        data_derivative = rho*projected(T_p+1:end,3);
        partial_mat{ind_r}(:,first_data+(1:T_d)) = ...
            sparse(T_p+(1:T_d),1:T_d,data_derivative,T,T_d);
        partial_mat{ind_r}(:,first_data+T_d+(1:T_d)) = ...
            sparse(T_p+(1:T_d),1:T_d,1j*data_derivative,T,T_d);
    end
end

FIM_mat_eta = sparse(eta_size,eta_size);
for ind_r1 = 1:RIS_K
    for ind_r2 = 1:RIS_K
        FIM_mat_eta = FIM_mat_eta + 2*input.noise_C_inv(ind_r1,ind_r2)* ...
            (partial_mat{ind_r1}'*partial_mat{ind_r2});
    end
end
Trans_mat = sparse(input.Trans_mat);
FIM_sparse = Trans_mat*real(FIM_mat_eta)*Trans_mat.' + sparse(input.FIM_para_pos_pri);
% Preserve the public dense FIM and covariance outputs. Small systems use
% the original dense inverse. For long frames, reuse one sparse LU and solve
% in column blocks: no dense identity or dense cubic factorization is needed.
FIM_para_pos = full(FIM_sparse);
npar = size(FIM_sparse,1);
if npar <= 512
    CRB = inv(FIM_para_pos);
else
    factor = decomposition(FIM_sparse,'lu');
    CRB = zeros(npar,npar);
    for first = 1:64:npar
        cols = first:min(first+63,npar);
        rhs = full(sparse(cols,1:numel(cols),1,npar,numel(cols)));
        CRB(:,cols) = factor\rhs;
    end
end
CRB_diag = diag(CRB);
pos_size = (3+2*RIS_K)*User_N;
bound = reshape(CRB_diag(1:pos_size),3+2*RIS_K,User_N);
obj1 = sum(bound(1:3,:),'all')/User_N;
obj2 = sum(CRB_diag(pos_size+1:end))/User_N/T_d;
end
