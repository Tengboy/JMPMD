function a = linespec(theta,N)
    a = exp(1j*[0:N-1].'*theta.');
end