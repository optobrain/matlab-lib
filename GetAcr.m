
%% not return abs(xcr)
% bsame: same length for convolution?  [defalut=false] https://app.asana.com/0/653197950698987/644340299733166

function ret = GetAcr(RR, maxlag, bsame)

if nargin < 3
    bsame = false;
end

	[nz,nx,ny,nt] = size(RR);
	if maxlag > nt
		error('maxlag > nt');
	end

	ret = zeros(nz,nx,ny,maxlag+1);
	parfor iy=1:ny
        RRy = RR(:,:,iy,:);
		for it=1:maxlag+1
            if bsame
                t1 = (1:nt-maxlag);  t2 = t1+it-1;
                M = sqrt(mean(abs(RRy(:,:,1,t1)).^2,4)) .* sqrt(mean(abs(RRy(:,:,1,t2)).^2,4));
                if min(M(:)) > 0
                    ret(:,:,iy,it) = mean( conj(RRy(:,:,1,t1)) .* RRy(:,:,1,t2) ,4) ./ M;
                else
                    ret(:,:,iy,it) = 0;
                end
            else
                t1 = (1:nt-it+1);  t2 = (it:nt);
                ret(:,:,iy,it) = mean( RRy(:,:,1,t2) .* conj(RRy(:,:,1,t1)) ,4) ./ mean( RRy(:,:,1,t1).*conj(RRy(:,:,1,t1)) ,4);
            end
		end
    end
   