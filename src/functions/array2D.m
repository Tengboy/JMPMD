function a = array2D(theta1,theta2,Nx,Ny)
    ax = exp(1j*[0:Nx-1].'*theta1.');
    ay = exp(1j*[0:Ny-1].'*theta2.');
    a = KR(ax,ay);
end