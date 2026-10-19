function [HH,Hmean] = vb_data_gram(omega,ax,ay,gainMean,gainCov,precision)
%VB_DATA_GRAM E[H'*precision*H] under the existing mean-field factors.
% Independent RIS gain vectors, independent per-user/axis angles; each
% gain vector retains its FULL user covariance. No independent-user gain
% approximation is made here. precision is fixed and may be non-diagonal.
% Shapes: omega=N-by-T-by-M; ax=Nx-by-K-by-M; gainMean=K-by-M;
% gainCov=K-by-K-by-M; HH=K-by-K-by-T; Hmean=M-by-K-by-T.
[K,M] = size(gainMean);
T = size(omega,2);
Hmean = zeros(M,K,T);
correction = zeros(K,K,T);
for m=1:M
    [psi,powerPsi] = vb_steering_moments(omega(:,:,m),ax(:,:,m),ay(:,:,m));
    gainSecond = conj(gainMean(:,m))*gainMean(:,m).' + gainCov(:,:,m).';
    for t=1:T
        h = gainMean(:,m).'.*psi(t,:);
        Hmean(m,:,t) = h;
        pp = psi(t,:)'*psi(t,:);
        pp(1:K+1:end) = powerPsi(t,:);
        % Across distinct RISs, E[conj(H_mi)*H_nj] factors into means.
        correction(:,:,t) = correction(:,:,t) + ...
            precision(m,m)*(gainSecond.*pp-h'*h);
    end
end
HH = zeros(K,K,T);
for t=1:T
    Q = Hmean(:,:,t)'*precision*Hmean(:,:,t)+correction(:,:,t);
    HH(:,:,t) = (Q+Q')/2;
end
end
