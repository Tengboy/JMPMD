function [theta,diagnostic]=pilot_initial_angles(cfg,observed)
% One-time, prior-local pilot initializer; no true channels/angles/data are read.
% Profile the unknown complex gains and the other users' pilot columns.
% Per-RIS marginal noise variances are used only by this initialization;
% the subsequent estimator retains its existing full noise covariance.
K=cfg.User_k;M=cfg.RIS_m;Tp=cfg.T_p;nx=cfg.RIS_size_x;ny=cfg.RIS_size_y;
assert(Tp>K && size(observed.X_p,1)==Tp,'Scaling:TooFewPilots','More pilot slots than users are required.');
theta=zeros(K,M,2);variance=theta;
for m=1:M
    for k=1:K
        offset=observed.User_pos_c(:,k)-observed.RIS_pos(:,m);
        distance=norm(offset);direction=offset/distance;
        for axis=1:2
            if axis==1,axisVector=observed.RIS_dir_x(:,m);else,axisVector=observed.RIS_dir_y(:,m);end
            theta(k,m,axis)=pi*axisVector.'*direction;
            gradient=pi*(axisVector-direction*(direction.'*axisVector))/distance;
            variance(k,m,axis)=max(cfg.pos_var*sum(gradient.^2),1e-12);
        end
    end
end
prior=theta;kappa=1./variance;sigma=sqrt(variance);
gram=observed.A_BS'*observed.A_BS;
y=((gram\observed.A_BS')*observed.Y(:,1:Tp)).';
covariance=observed.effective_noise_power*(gram\eye(M));nu=real(diag(covariance));
assert(all(isfinite(nu)&nu>0));
initial=zeros(1,M);final=initial;
for m=1:M
    W=observed.Omega(:,1:Tp,m).';Xp=observed.X_p;
    D=Xp.*(W*array2D(theta(:,m,1),theta(:,m,2),nx,ny));
    initial(m)=objective(D,y(:,m),nu(m),theta(:,m,:),prior(:,m,:),kappa(:,m,:));
    for sweep=1:3
        for k=1:K
            other=D(:,[1:k-1,k+1:K]);[U,S,~]=svd(other,'econ');singular=diag(S);
            Q=U(:,singular>max(size(other))*eps(max(singular)));
            residual=y(:,m)-Q*(Q'*y(:,m));
            for level=1:3
                center=reshape(theta(k,m,:),1,2);origin=reshape(prior(k,m,:),1,2);
                sd=reshape(sigma(k,m,:),1,2);concentration=reshape(kappa(k,m,:),1,2);
                if level==1
                    gx=origin(1)+(-4:4)*(.75*sd(1));gy=origin(2)+(-4:4)*(.75*sd(2));
                else
                    radius=3/4^(level-1);
                    gx=center(1)+(-2:2)*(radius*sd(1)/2);gy=center(2)+(-2:2)*(radius*sd(2)/2);
                end
                gx=unique(min(pi,max(-pi,min(origin(1)+3*sd(1),max(origin(1)-3*sd(1),gx)))));
                gy=unique(min(pi,max(-pi,min(origin(2)+3*sd(2),max(origin(2)-3*sd(2),gy)))));
                [xx,yy]=ndgrid(gx,gy);candidates=[xx(:),yy(:);center;origin];
                responses=W*array2D(candidates(:,1),candidates(:,2),nx,ny);
                columns=Xp(:,k).*responses;projected=columns-Q*(Q'*columns);
                energy=real(sum(abs(projected).^2,1));
                penalty=sum(concentration.*(1-cos(candidates-origin)),2).';
                score=abs(residual'*projected).^2./max(nu(m)*energy,realmin)-penalty;
                [~,best]=max(score);theta(k,m,:)=candidates(best,:);D(:,k)=columns(:,best);
            end
        end
    end
    final(m)=objective(D,y(:,m),nu(m),theta(:,m,:),prior(:,m,:),kappa(:,m,:));
end
diagnostic=struct('prior_theta',prior,'prior_variance',variance,'initial_objective',initial, ...
    'final_objective',final,'sweeps',3,'search_sigma',3,'used_slots',1:Tp);
end
function value=objective(D,y,nu,theta,prior,kappa)
residual=y-D*(D\y);
value=real(residual'*residual)/nu+sum(kappa.*(1-cos(theta-prior)),'all');
end
