%% Get NORMALIZED Gaussian convolution kernel
% sig = exp( -x^2 / 2 / sig^2 )
% As default, GaussKernel([5 5 5]) returns a kernel of sig=1 (FWHM=2.35)

function res = GaussKernel(arySize, sig)

if nargin < 2
	sig = (arySize-1)/4 +eps;		% boundary will have exp(-2^2/2*1^2) = 0.14
end

if size(arySize) ~= size(sig)
	error('arySize and sig have different sizes.');
end
	
	cnv = zeros(arySize);
	cnt = (arySize+1)/2;		% center

	switch numel(arySize)
		case 1
			for i1=1:arySize(1)
				cnv(i1) = exp( -(i1-cnt)^2 / 2 / sig^2 );
			end
		case 2
			for i1=1:arySize(1)
				for i2=1:arySize(2)
					cnv(i1,i2) = exp( -norm(([i1 i2]-cnt)./sig)^2 / 2 );
				end
			end
		case 3
			for i1=1:arySize(1)
				for i2=1:arySize(2)
					for i3=1:arySize(3)
						cnv(i1,i2,i3) = exp( -norm(([i1 i2 i3]-cnt)./sig)^2 / 2 );
					end
				end
			end
		case 4
			for i1=1:arySize(1)
				for i2=1:arySize(2)
					for i3=1:arySize(3)
						for i4=1:arySize(4)
							cnv(i1,i2,i3,i4) = exp( -norm(([i1 i2 i3 i4]-cnt)./sig)^2 / 2 );
						end
					end
				end
			end
		case 5
			for i1=1:arySize(1)
				for i2=1:arySize(2)
					for i3=1:arySize(3)
						for i4=1:arySize(4)
							for i5=1:arySize(5)
								cnv(i1,i2,i3,i4,i5) = exp( -norm(([i1 i2 i3 i4 i5]-cnt)./sig)^2 / 2 );
							end
						end
					end
				end
			end
	end

	res = cnv / sum(cnv(:));

    