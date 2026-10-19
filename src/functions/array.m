function a_r=array(theta,N)
    a_r = exp(1j*[0:N-1].'*theta.');
end