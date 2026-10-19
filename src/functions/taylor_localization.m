function [ue_pos,ue_pos_C,flag,diagnostics] = taylor_localization(measure,para,cov,initpos)
%TAYLOR_LOCALIZATION Weighted Taylor position update, paper (29)-(30).
% measure and cov are RIS-count by 2, ordered [theta_x,theta_y].
% Angles use spatial-frequency radians: theta = pi * direction cosine.
% The third output retains its historical meaning: iteration count, NOT a
% convergence flag. The caller must decide how to handle a failed estimate.
measurementVariance = reshape(cov.',[],1);
ue_pos = initpos;
delta = ones(3,1);
flag = 0;

% Keep the original tolerance, step cutoff and maximum iteration count.
while norm(delta)>=1e-8 && norm(delta)<100 && flag<30
    [A,residual] = linearize_geometry(ue_pos,measure,para);
    weightedA = A./measurementVariance;
    information = A.'*weightedA;
    delta = information\(weightedA.'*residual);
    ue_pos = ue_pos+delta;
    flag = flag+1;
end

% Evaluate covariance at the returned position, including at an early cutoff.
A = linearize_geometry(ue_pos,measure,para);
information = A.'*(A./measurementVariance);
ue_pos_C = information\eye(3);
% Optional diagnostics never replace, regularize, or discard the estimate.
if nargout >= 4
    diagnostics.iterations = flag;
    diagnostics.last_step_norm = norm(delta);
    diagnostics.information_rcond = rcond(information);
    finiteResult = all(isfinite(ue_pos(:))) && all(isfinite(ue_pos_C(:))) ...
        && isfinite(diagnostics.last_step_norm);
    % reason_code: 0 converged; 1 nonfinite; 2 step cutoff; 3 iteration cap;
    % 4 ill-conditioned final information; 5 invalid measurement variance.
    if ~finiteResult
        code = 1;
    elseif any(isnan(measurementVariance)) || any(measurementVariance<=0)
        code = 5;
    elseif diagnostics.information_rcond < 1e-12
        code = 4;
    elseif diagnostics.last_step_norm >= 100
        code = 2;
    elseif diagnostics.last_step_norm >= 1e-8
        code = 3;
    else
        code = 0;
    end
    diagnostics.reason_code = code;
    diagnostics.converged = code==0;
end
end

function [A,residual] = linearize_geometry(position,measure,para)
numberOfAnchors = size(para.RIS_pos,2);
A = zeros(2*numberOfAnchors,3);
residual = zeros(2*numberOfAnchors,1);
for anchor = 1:numberOfAnchors
    offset = position-para.RIS_pos(:,anchor);
    distance = norm(offset);
    direction = offset/distance;
    axes = [para.RIS_dir_x(:,anchor),para.RIS_dir_y(:,anchor)];
    rows = 2*anchor-1:2*anchor;
    projection = axes.'*direction;
    % d(pi*e^T*u)/dp = pi/r * (e^T - (e^T*u)*u^T).
    A(rows,:) = (pi/distance)*(axes.'-projection*direction.');
    residual(rows) = measure(anchor,:).'-pi*projection;
end
end
