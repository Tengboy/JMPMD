function lagSum = ris_lag_sum(R, nx, ny)
%RIS_LAG_SUM Sum UPA correlation entries with the same spatial lag.
% KR(ax,ay) stores y fastest. The output therefore has y lags in rows,
% x lags in columns, and one page per RIS. Supports rectangular UPAs.
persistent cachedShape lagIndex
if isempty(cachedShape) || ~isequal(cachedShape, [nx, ny])
    element = (0:nx*ny-1).';
    x = floor(element/ny);
    y = mod(element,ny);
    index = (y.'-y+ny) + (2*ny-1)*(x.'-x+nx-1);
    % Match the original i-outer/j-inner accumulation order.
    lagIndex = reshape(index.', [], 1);
    cachedShape = [nx, ny];
end
lagSum = zeros(2*ny-1, 2*nx-1, size(R,3), 'like', R);
for ris = 1:size(R,3)
    entries = reshape(R(:,:,ris).', [], 1);
    sums = accumarray(lagIndex, entries, [(2*ny-1)*(2*nx-1), 1], @sum);
    lagSum(:,:,ris) = reshape(sums, 2*ny-1, 2*nx-1);
end
end
