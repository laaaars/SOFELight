function R = reshapeTop(vec, slots, varargin)
% creates matrix with column i containing the next slots(i) entries from vec
% Columnwise the matrix is filled from top to bottom
% remaining entries are filled with zero

if sum(slots) ~= numel(vec)
    error('! Incompatible data !');
end
if nargin > 2
    cols = varargin{1};
    if cols < max(slots)
        error('! Dimension mismatch!');
    end
else
    cols = max(slots);
end
rows = numel(slots);
R = repmat((1:cols)',1,rows); 
R = bsxfun(@le,R,slots(:)');
[ii,jj] = find(R>0);
R = full(sparse(ii, jj, vec, cols, rows));