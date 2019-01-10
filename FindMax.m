%% Find maximum for multi-dimensional (2 - 5) array
%% 03/25/10, Jonghwan

% ary = array, ndims(ary) = 2 ~ 5

function [val, idx] = FindMax(ary, bSigned)
if nargin < 2
	bSigned = 0;
end

	ary = squeeze(ary);
	if (ndims(ary) > 6)
		error('Dimension of ary is larger than 5.');
	end

	val = 0;
	idx = zeros(1,ndims(ary));
	if bSigned == 0
		switch (ndims(ary))
			case (1)
				[val,idx] = max(ary);
			case (2)
				[m,i1] = max(ary);
				[val,i2] = max(m);
				idx(2) = i2;
				idx(1) = i1(1,idx(2));
			case (3)
				[m,i1] = max(ary);
				[m,i2] = max(m);
				[val,i3] = max(m);
				idx(3) = i3;
				idx(2) = i2(1,1,idx(3));
				idx(1) = i1(1,idx(2),idx(3));
				
			case (4)
				[m,i1] = max(ary);
				[m,i2] = max(m);
				[m,i3] = max(m);
				[val,i4] = max(m);
				idx(4) = i4;
				idx(3) = i3(1,1,1,idx(4));
				idx(2) = i2(1,1,idx(3),idx(4));
				idx(1) = i1(1,idx(2),idx(3),idx(4));
				
			case (5)
				[m,i1] = max(ary);
				[m,i2] = max(m);
				[m,i3] = max(m);
				[m,i4] = max(m);
				[val,i5] = max(m);
				idx(5) = i5;
				idx(4) = i4(1,1,1,1,idx(5));
				idx(3) = i3(1,1,1,idx(4),idx(5));
				idx(2) = i2(1,1,idx(3),idx(4),idx(5));
				idx(1) = i1(1,idx(2),idx(3),idx(4),idx(5));
				
			case (6)
				[m,i1] = max(ary);
				[m,i2] = max(m);
				[m,i3] = max(m);
				[m,i4] = max(m);
				[m,i5] = max(m);
				[val,i6] = max(m);
				idx(6) = i6;
				idx(5) = i5(1,1,1,1,1,idx(6));
				idx(4) = i4(1,1,1,1,idx(5),idx(6));
				idx(3) = i3(1,1,1,idx(4),idx(5),idx(6));
				idx(2) = i2(1,1,idx(3),idx(4),idx(5),idx(6));
				idx(1) = i1(1,idx(2),idx(3),idx(4),idx(5),idx(6));
		end

	else
		switch (ndims(ary))
			case (1)
				[val,idx] = FindSignedMax(ary);
			case (2)
				[m,i1] = FindSignedMax(ary);
				[val,i2] = FindSignedMax(m);
				idx(2) = i2;
				idx(1) = i1(1,idx(2));
			case (3)
				[m,i1] = FindSignedMax(ary);
				[m,i2] = FindSignedMax(m);
				[val,i3] = FindSignedMax(m);
				idx(3) = i3;
				idx(2) = i2(1,1,idx(3));
				idx(1) = i1(1,idx(2),idx(3));
				
			case (4)
				[m,i1] = FindSignedMax(ary);
				[m,i2] = FindSignedMax(m);
				[m,i3] = FindSignedMax(m);
				[val,i4] = FindSignedMax(m);
				idx(4) = i4;
				idx(3) = i3(1,1,1,idx(4));
				idx(2) = i2(1,1,idx(3),idx(4));
				idx(1) = i1(1,idx(2),idx(3),idx(4));
				
			case (5)
				[m,i1] = FindSignedMax(ary);
				[m,i2] = FindSignedMax(m);
				[m,i3] = FindSignedMax(m);
				[m,i4] = FindSignedMax(m);
				[val,i5] = FindSignedMax(m);
				idx(5) = i5;
				idx(4) = i4(1,1,1,1,idx(5));
				idx(3) = i3(1,1,1,idx(4),idx(5));
				idx(2) = i2(1,1,idx(3),idx(4),idx(5));
				idx(1) = i1(1,idx(2),idx(3),idx(4),idx(5));
				
			case (6)
				[m,i1] = FindSignedMax(ary);
				[m,i2] = FindSignedMax(m);
				[m,i3] = FindSignedMax(m);
				[m,i4] = FindSignedMax(m);
				[m,i5] = FindSignedMax(m);
				[val,i6] = FindSignedMax(m);
				idx(6) = i6;
				idx(5) = i5(1,1,1,1,1,idx(6));
				idx(4) = i4(1,1,1,1,idx(5),idx(6));
				idx(3) = i3(1,1,1,idx(4),idx(5),idx(6));
				idx(2) = i2(1,1,idx(3),idx(4),idx(5),idx(6));
				idx(1) = i1(1,idx(2),idx(3),idx(4),idx(5),idx(6));
		end

	end


