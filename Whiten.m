
% Y(ch,t)
% bnorm = normalize by std?

function ret = Whiten(Y,bnorm)

if nargin < 2
    bnorm = 1;
end

    [nch nt] = size(Y);
    
    if bnorm
        ret = ( Y - repmat(mean(Y,2),[1 nt]) ) ./ repmat(std(Y,1,2),[1 nt]);
    else
        ret = Y - repmat(mean(Y,2),[1 nt]);
    end
    
    