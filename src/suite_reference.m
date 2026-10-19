function ref=suite_reference(cfg,scene)
% Existing realization-wise inverse information bounds, then MC averaging.
% Full-covariance joint MAP is the single published QPSK ideal-CSI curve.
% The diagonal metric is retained only as a diagnostic, not plotted.
% The raw FIM is poorly scaled because gains are ~1e-10. Store raw/scaled
% condition and inverse residual for EVERY reference; avoid millions of
% repeated scale warnings. This does not change the existing inverse.
near=warning('query','MATLAB:nearlySingularMatrix');
restore=onCleanup(@()warning(near)); %#ok<NASGU>
warning('off','MATLAB:nearlySingularMatrix');
if cfg.kind=='Q', ref=reference_q(cfg,scene); else, ref=reference_g(cfg,scene); end
ref.csi_ser=NaN; ref.csi_full_ser=NaN;
if cfg.kind=='Q'
    states=scene.constellation(Constell_Mat(cfg.User_k,4));
    projection=(scene.A_BS'*scene.A_BS)\scene.A_BS';
    y=projection*scene.Y; cov=scene.noise_power*inv(scene.A_BS'*scene.A_BS);
    precisions=cat(3,diag(1./diag(cov)),inv(cov));
    errors=zeros(1,2);
    for mode=1:2
        for t=1:scene.T_d
            H=scene.H(:,:,t); z=y(:,cfg.T_p+t); precision=precisions(:,:,mode);
            residual=z-H*states.';
            costs=real(sum(conj(residual).*(precision*residual),1));
            [~,idx]=min(costs);
            errors(mode)=errors(mode)+nnz(states(idx,:)~=scene.X_d(t,:));
        end
    end
    ref.csi_ser=errors(1)/(cfg.User_k*scene.T_d);
    ref.csi_full_ser=errors(2)/(cfg.User_k*scene.T_d);
end
end
