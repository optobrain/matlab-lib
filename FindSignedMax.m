
function [m, im] = FindSignedMax(A, dim)
if nargin < 2
	dim = 1;
end

%{
	for id=1:ndims(A)
		if size(A,id) == 1
			if dim == id
				dim = id + 1;
			end
		else
			break;
		end
	end
%}
	if min(A(:)) < 0
		m0 = max(max(A,[],dim),0);
		[mA im] = max(abs(A),[],dim);
		m = mA .* (1-2*ceil((mA-m0)./max(mA,eps)));
	else
		[m im] = max(A,[],dim);
	end
