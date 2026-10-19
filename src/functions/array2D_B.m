function a = array2D_B(theta1,theta2,kappa1,kappa2,Nx,Ny)
    ax = exp(1j*[0:Nx-1].'*theta1.').*( besseli([0:Nx-1].',kappa1,1)/besseli(0,kappa1,1) );
    ay = exp(1j*[0:Ny-1].'*theta2.').*( besseli([0:Ny-1].',kappa2,1)/besseli(0,kappa2,1) );
    a = KR(ax,ay);
end