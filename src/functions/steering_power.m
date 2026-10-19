function value=steering_power(omega,ax,ay)
% E|omega^T a|^2 for separable UPA moments, without an N-by-N matrix.
% vec(By*V*Bx.') = kron(Bx,By)*vec(V); y index varies fastest.
nx=size(ax,1); ny=size(ay,1); T=size(omega,2);
value=zeros(T,size(ax,2));
for k=1:size(ax,2)
    Bx=toeplitz(ax(:,k),ax(:,k)'); By=toeplitz(ay(:,k),ay(:,k)');
    v=reshape(By*reshape(conj(omega),ny,[]),ny,nx,T);
    v=Bx*reshape(permute(v,[2 1 3]),nx,[]);
    v=reshape(permute(reshape(v,nx,ny,T),[2 1 3]),nx*ny,T);
    value(:,k)=real(sum(omega.*v,1)).';
end
end
