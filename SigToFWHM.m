% sigma of canonical Guassian shape to FWHM
% exp( - x^2 /2 /sig^2 )

function ret = SigToFWHM(sig)

	ret = 2*sqrt(2*log(2))*sig;


	
		
