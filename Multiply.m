%% Multiply n-dim arrays
%% 05/30/10, Jonghwan

function ret = Multiply(A, B, bCommon)
if nargin < 3
	bCommon = 1;
end
	% bCommon = 1 : A and B has common index
	%	eg. A = [2 3 5], B = [5 6] will return [2 3 6]
	% bCommon = 0 : eg. A = [2 3 5], B = [5 6] will return [2 3 5 5 6]

	szA = size(A);
	szB = size(B);

	if bCommon == 1
		if (szA(end) ~= szB(1))
			error('Two arrays should have common dimension when bCommon=1.');
		end
		szret = [szA(1:end-1) szB(2:end)];
		ret = zeros([numel(A)/szA(end) numel(B)/szB(1)]);

		sqA = reshape(A, [numel(A)/szA(end) szA(end)]);
		sqB = reshape(B, [szB(1) numel(B)/szB(1)]);

		for (j=1:szA(end))
			ret = ret + sqA(:,j) * sqB(j,:);
		end

		ret = reshape(ret,szret);

	else
		szret = [szA szB];
		if ndims(szB) == 2 && szB(1) == 1
			szret = [szA szB(2)];
		end
		sqA = reshape(A, [numel(A) 1]);
		sqB = reshape(B, [1 numel(B)]);

		ret = sqA * sqB;
		ret = reshape(ret,szret);
	end
