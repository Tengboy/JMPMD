function [stableCount,relativeChange] = iteration_status(current,previous,stableCount,tolerance)
%ITERATION_STATUS Compare estimates only; no simulation truth is an input.
% Every supplied state field must be stable for two successive iterations
% (the caller decides the required count). Scale by both current/old norms;
% realmin, rather than eps, preserves sensitivity to small physical gains.
relativeChange = 0;
names = fieldnames(current);
for n=1:numel(names)
    a = current.(names{n}); b = previous.(names{n});
    if any(~isfinite(a(:))) || any(~isfinite(b(:)))
        relativeChange = Inf;
        break;
    end
    scale = max([norm(a(:)),norm(b(:)),realmin]);
    relativeChange = max(relativeChange,norm(a(:)-b(:))/scale);
end
if relativeChange < tolerance
    stableCount = stableCount+1;
else
    stableCount = 0;
end
end
