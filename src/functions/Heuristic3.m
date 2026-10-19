function [a, mu, kappa] = Heuristic3(eta, m, prior_mu, prior_ka, ini)
% Approximate the angular likelihood and VM prior by one VM distribution.
eta = eta(:);
m = m(:);
if numel(eta) ~= numel(m)
    % Existing experiment callers retain full eta but pass only nonzero lags.
    % Restore the complete lag vector, including trailing zero coefficients.
    if isequal(m, find(abs(eta)>0)-1)
        m = (0:numel(eta)-1)';
    else
        error('Heuristic3:SizeMismatch','eta and m must have matching lengths.');
    end
end
if isempty(m) || any(~isfinite(eta)) || any(~isfinite(m)) || ...
        any(m<0 | m~=fix(m)) || ~isscalar(prior_ka) || ...
        ~isfinite(prior_ka) || prior_ka<0 || ~isfinite(prior_mu) || ~isfinite(ini)
    error('Heuristic3:InvalidInput','Inputs must be finite, with nonnegative integer lags and nonnegative prior concentration.');
end
n = (0:max(m))';
eta_q0 = prior_ka*exp(1j*prior_mu);
if ~any(eta(m>0))
    % A zero or constant likelihood leaves the prior exactly unchanged.
    kappa = prior_ka;
    mu = prior_mu;
    if kappa==0
        mu=0;
    end
    a = vm_moments(n,mu,kappa);
    return;
end

mu = ini;
m2 = m.^2;
gradient3 = @(x)real(1j*eta_q0' * exp(1j*x) + 1j*eta' * (m.*exp(1i*m*x)));
curvature = @(x)-real(eta_q0' * exp(1i*x) + eta' * (m2.*exp(1i*m*x)));
objv = @(x)real(eta_q0' * exp(1i*x) + eta' * exp(1i*m*x));
% Follow the local mode selected by the current iterate. A full-circle grid
% can jump to an unrelated lobe when the local curvature is nonnegative.
% Limit phase motion at the largest lag, then accept only uphill steps.
max_step = pi/(4*max(1,max(m)));
gradient_scale = max(1,abs(eta_q0)+sum(abs(eta).*m));
for count=1:200
    d=gradient3(mu);
    d2=curvature(mu);
    if d2<0 && abs(d)<=1e-12*gradient_scale
        break;
    end
    if d2<0
        direction=max(-max_step,min(max_step,-d/d2));
    else
        direction=sign(d)*max_step;
        if abs(d)<=1e-12*gradient_scale
            % Escape an exactly stationary minimum using local probes only.
            probe=min(max_step,1e-3);
            if objv(mu+probe)>=objv(mu-probe), direction=probe;
            else, direction=-probe; end
        end
    end
    old_value=objv(mu); step=1; accepted=false;
    for backtrack=1:50
        mu_n=mu+step*direction;
        if objv(mu_n)>=old_value
            accepted=true;
            break;
        end
        step=step/2;
    end
    if ~accepted, break; end
    delta=abs(mu_n-mu);
    mu=mu_n;
    if delta<1e-10, break; end
end
d2=curvature(mu);
if ~isfinite(d2) || d2>=0
    % A flat peak has no finite positive curvature-based concentration.
    kappa=0;
else
    moment=exp(0.5/d2);
    if moment==1
        % Avoid Ainv(1)=Inf when a sharply concentrated moment rounds to 1.
        kappa=-d2;
    else
        kappa=Ainv(moment);
    end
end
a=vm_moments(n,mu,kappa);
end

function a=vm_moments(n,mu,kappa)
if kappa==0
    a=zeros(size(n));
    a(n==0)=1;
else
    a=exp(1i*n*mu).*(besseli(n,kappa,1)/besseli(0,kappa,1));
end
end
