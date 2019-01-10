% Num2Str(2,3) = '002'

function ret = NumToStr(num, nstr)


	if nargin < 2
		ret = num2str(num);
	else
		ret = num2str(num);
		if nstr > length(ret)
			for ii=1:nstr-length(ret)
				ret = ['0' ret];
			end
		end
	end

	
		
