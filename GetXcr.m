%% Get the COMPLEX-valued cross-correlation of time-series volume data RR relative to R0

% RR [nz nx ny nt]
% R0 [nz nx ny] : reference volume data
% bmag : xcr of magnitude of RR : RR = RR - mean(RR,t)


function ret = GetXcr(RR, R0, bmag)

	% argin
		if nargin < 3
			bmag = 0;
		end
		if nargin < 2
			R0 = RR(:,:,:,1);
		end

	% xcr
		if bmag == 1
			RR = abs(RR);
			RR = RR - repmat(mean(mean(mean(RR,1),2),3),[size(RR,1) size(RR,2) size(RR,3) 1]);
			R0 = abs(R0);
			R0 = R0 - mean(mean(mean(R0,1),2),3);
			ret = sum(sum(sum( RR .* repmat(R0,[1 1 1 size(RR,4)]) ,1),2),3) ./ sqrt(sum(sum(sum( RR.^2 ,1),2),3)) ./ repmat( sqrt(sum(sum(sum( R0.^2 ,1),2),3)) ,[1 1 1 size(RR,4)]);
		else
			ret = sum(sum(sum( RR .* repmat(conj(R0),[1 1 1 size(RR,4)]) ,1),2),3) ./ sqrt(sum(sum(sum( abs(RR).^2 ,1),2),3)) ./ repmat( sqrt(sum(sum(sum( abs(R0).^2 ,1),2),3)) ,[1 1 1 size(RR,4)]);
		end
	
