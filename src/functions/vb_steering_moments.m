function [meanPsi,powerPsi] = vb_steering_moments(omega,ax,ay)
%VB_STEERING_MOMENTS E[w^T a] and E[|w^T a|^2] under q(theta_x)q(theta_y).
% omega: (Nx*Ny)-by-T; ax/ay: E[exp(1j*n*theta)] with n=0,...,N-1.
% Users are independent under the current variational factorization.
% Explicit first column/row avoids ambiguity in complex toeplitz(v).
meanPsi = omega.'*KR(ax,ay);
powerPsi = steering_power(omega,ax,ay);
end
