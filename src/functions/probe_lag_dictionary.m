function dictionary=probe_lag_dictionary(omega,nx,ny)
% Column t contains sum conj(w_i,t)*w_j,t for each UPA lag j-i.
% Zero padding gives linear (not circular) correlation. y varies fastest.
grid=reshape(omega,ny,nx,[]);
frequency=fft2(grid,2*ny-1,2*nx-1);
correlation=fftshift(fftshift(ifft2(conj(frequency).*frequency),1),2);
dictionary=reshape(correlation,(2*ny-1)*(2*nx-1),[]);
end
