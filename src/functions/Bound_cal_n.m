% Position information bound with every transmitted symbol known.
function [obj,FIM_para_pos,CRB] = Bound_cal_n(input)
User_N = size(input.UE_location_r,2);
RIS_K = size(input.RIS_para.RIS_pos,2);
Nx = input.RIS_para.RIS_size_x;
Ny = input.RIS_para.RIS_size_y;
block_size = 4*User_N;
eta_size = block_size*RIS_K;
coef = zeros(1,eta_size);
projected = cell(RIS_K,1);

% Project each RIS/user derivative once, rather than rebuilding it for every
% pair of RISs and users in the Fisher-matrix cross products.
for ind_r = 1:RIS_K
    A = zeros(Nx*Ny,block_size);
    for ind_u = 1:User_N
        first = (ind_u-1)*4;
        global_first = (ind_r-1)*block_size + first;
        rho = input.rho(ind_u,ind_r);
        coef(global_first+(1:4)) = [rho,rho,rho/abs(rho),1j*rho];
        ax = linespec(input.th_RIS_x(ind_u,ind_r),Nx);
        ay = linespec(input.th_RIS_y(ind_u,ind_r),Ny);
        a = KR(ax,ay);
        A(:,first+(1:4)) = [KR((1j*(0:Nx-1).').*ax,ay), ...
            KR(ax,(1j*(0:Ny-1).').*ay),a,a];
    end
    projected{ind_r} = input.X(:,repelem(1:User_N,4)).* ...
        (input.Omega1(:,:,ind_r).'*A);
end
Lambda = zeros(eta_size,eta_size);
for ind_r1 = 1:RIS_K
    rows = (ind_r1-1)*block_size+(1:block_size);
    for ind_r2 = 1:RIS_K
        cols = (ind_r2-1)*block_size+(1:block_size);
        Lambda(rows,cols) = 2*input.noise_C_inv(ind_r1,ind_r2)* ...
            (projected{ind_r1}'*projected{ind_r2});
    end
end
Trans_mat = input.Trans_mat;
FIM_para_pos = full(Trans_mat*real(Lambda.*(coef'*coef))*Trans_mat.' + input.FIM_para_pos_pri);
CRB = inv(FIM_para_pos);
bound = reshape(diag(CRB),3+2*RIS_K,User_N);
obj = sum(bound(1:3,:),'all')/User_N;
end
