
function [res, clim] = Rescale(img, toClim, fromClim, bforce)
if nargin < 4
    bforce = 1;
end

if nargin < 3
	fromClim = zeros(1,2);
	fromClim(1) = min(img(:));
	fromClim(2) = max(img(:));
else
    sz = size(img);
	img = img(:);
    if bforce == 1
        img(img<fromClim(1)) = fromClim(1);
        img(img>fromClim(2)) = fromClim(2);
    end
    img = reshape(img,sz);
end
	clim = [min(img(:)) max(img(:))];
	res = (toClim(2)-toClim(1))/(fromClim(2)-fromClim(1)) * (img - fromClim(1)) + toClim(1);
