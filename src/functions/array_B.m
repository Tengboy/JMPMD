function a_r=array_B(theta,kappa,N)
    a_r = exp(1j*[0:N-1].'*theta.').*( besseli([0:N-1].',kappa,1)/besseli(0,kappa,1) );
end