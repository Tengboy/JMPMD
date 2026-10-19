function test_strict_moments(sourceDir)
% Independent finite-distribution sums and circular quadrature, not copies
% of the implementation. Covariances use E[(z-mu)*(z-mu)'] throughout.
addpath(fullfile(sourceDir,'functions'),'-begin');
assert(exist('vb_steering_moments','file')==2, 'Strict moment implementation missing.');
rng(291);
nx=3; ny=2; T=4; K=2; M=2;
omega=randn(nx*ny,T)+1j*randn(nx*ny,T);
mu=[0.7,-0.9]; kappa=[1.2,3.1];
ax=array_B(mu.',kappa(1),nx); % array_B takes angle columns.
muY=[-0.4,1.1]; ay=array_B(muY.',kappa(2),ny);
[pm,p2]=vb_steering_moments(omega,ax,ay);
grid=-pi+(0:1023)'*2*pi/1024;
for k=1:K
    px=exp(kappa(1)*cos(grid-mu(k))); px=px/sum(px);
    py=exp(kappa(2)*cos(grid-muY(k))); py=py/sum(py);
    sx=exp(1j*(0:nx-1)'*grid.'); sy=exp(1j*(0:ny-1)'*grid.');
    % Integrate over x numerically; integrate y using explicit quadrature
    % covariance. Neither quadrature uses Bessel or Toeplitz constructions.
    my=sy*py; By=(sy.*py.')*sy';
    meanExact=zeros(T,1); powerExact=zeros(T,1);
    for g=1:numel(grid)
        W=reshape(omega,ny,nx,T);
        projectedY=zeros(ny,T);
        for t=1:T, projectedY(:,t)=W(:,:,t)*sx(:,g); end
        meanExact=meanExact+px(g)*(projectedY.'*my);
        powerExact=powerExact+px(g)*real(sum(projectedY.*(By*conj(projectedY)),1)).';
    end
    assert(norm(pm(:,k)-meanExact)<1e-10);
    assert(norm(p2(:,k)-powerExact)<1e-10);
end
fprintf('PASS VM array moments against circular quadrature (rectangular 3x2).\n');

% Two correlated users in each discrete factor; independent RIS factors.
xs=[1+1j,2-1j; -2+0.3j,0.5+2j; 0.2-1j,-1-0.7j]; px=[.2;.5;.3];
as=[0.7+1j,-0.4+0.2j; -1+0.4j,2+0.1j]; pa=[.35;.65];
[xm,xc]=discrete_posterior_moments(xs,px);
[am,ac]=discrete_posterior_moments(as,pa);
% Steering factor has two independent states per user, exact finite moments.
s1=[1+2j,-0.6+0.1j]; s2=[-.3+.7j,1.5-1j]; ps=[.4,.6];
sm=[s1*ps.',s2*ps.']; sp=[abs(s1).^2*ps.',abs(s2).^2*ps.'];
Jexact=zeros(K); Qexact=zeros(K);
for ix=1:3
    for i=1:2
        for j=1:2
            a=xs(ix,:).*[s1(i),s2(j)];
            Jexact=Jexact+px(ix)*ps(i)*ps(j)*(a'*a);
        end
    end
end
J=vb_gain_gram(xm,reshape(xc,K,K,1),sm,sp);
assert(norm(J-Jexact,'fro')<1e-11);
fprintf('PASS gain Gram against joint finite-distribution enumeration.\n');

% Test full data Gram by enumerating independent gain vectors and independent
% angle factors for two RISs. Use one element per axis with a deterministic
% steering case separately below; the generic projected moments are tested
% through the data helper's second output and explicit angle grid integration.
Om=cat(3,omega,omega.*exp(.4j));
Ax=repmat(ax,1,1,M); Ay=repmat(ay,1,1,M);
gain=[am.',(am*[1,.2j;0,1]).'];
Cgain=cat(3,ac,[1,.2j;0,1].'*ac*conj([1,.2j;0,1]));
P=[1,.2+.1j;.2-.1j,1.4];
[Q,Hbar]=vb_data_gram(Om,Ax,Ay,gain,Cgain,P);
% Integrate channel Gram using high-accuracy steering quadrature moments,
% derived above, then sum over the finite gain distribution explicitly.
gainStates=cat(3,as,as*[1,.2j;0,1]);
for t=1:T
    Qexact=zeros(K);
    for ia=1:2
        for ib=1:2
            gainState=[gainStates(ia,:,1);gainStates(ib,:,2)];
            means=[pm(t,:);exp(.4j)*pm(t,:)];
            H=gainState.*means;
            q=H'*P*H;
            for m=1:M
                q=q+P(m,m)*diag(abs(gainState(m,:)).^2.*(p2(t,:)-abs(means(m,:)).^2));
            end
            Qexact=Qexact+pa(ia)*pa(ib)*q;
        end
    end
    assert(norm(Q(:,:,t)-Qexact,'fro')<1e-10);
    assert(norm(Q(:,:,t)-Q(:,:,t)','fro')<1e-11);
    assert(min(eig(Q(:,:,t)))>-1e-10);
end
fprintf('PASS data Gram with complex off-diagonal gain covariance and full precision.\n');

AxPoint=repmat(exp(1j*(0:nx-1)'*mu),1,1,M);
AyPoint=repmat(exp(1j*(0:ny-1)'*[-.4,1.1]),1,1,M);
[Qpoint,Hpoint]=vb_data_gram(Om,AxPoint,AyPoint,gain,zeros(K,K,M),P);
for t=1:T
    assert(norm(Qpoint(:,:,t)-Hpoint(:,:,t)'*P*Hpoint(:,:,t),'fro')<1e-10);
end
assert(all(isfinite(Hbar(:))));
fprintf('PASS deterministic-angle/zero-covariance limit.\n');
end
