function J = vb_gain_gram(xMean,xCov,meanPsi,powerPsi)
%VB_GAIN_GRAM E[A'*A], A(t,k)=x(t,k)*psi(t,k), for one RIS.
% xMean/meanPsi/powerPsi: T-by-K; xCov: K-by-K-by-T, with the
% standard complex convention Cov(i,j)=E[(x_i-mu_i)*conj(x_j-mu_j)].
% Therefore E[conj(x_i)*x_j] uses Cov(j,i), NOT Cov(i,j).
K = size(xMean,2);
J = zeros(K);
for t=1:size(xMean,1)
    xx = xMean(t,:)'*xMean(t,:) + xCov(:,:,t).';
    pp = meanPsi(t,:)'*meanPsi(t,:);
    pp(1:K+1:end) = powerPsi(t,:);
    J = J + xx.*pp;
end
J = (J+J')/2; % Remove roundoff only; do not clip eigenvalues.
end
