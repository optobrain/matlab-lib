%% plot image sequence
%% imgs should be cells containing images of the same size
%% alps = alpha datas

function ret = PlotImageSeq(imgs, clim, bgpix, bgclr, alps)

if nargin < 3
	bgpix = 1;
end
if nargin < 4
	bgclr = 'k';
end
if nargin < 5
	alps = imgs;
end

	[nimg1,nimg2] = size(imgs);
	[npix1,npix2] = size(imgs{1,1});

	img = zeros(nimg1*(npix1+bgpix)+bgpix,nimg2*(npix2+bgpix)+bgpix);
	alp = img;
	for iimg1=1:nimg1
		for iimg2=1:nimg2
%           [bgpix+(iimg1-1)*(npix1+bgpix)+1 iimg1*(npix1+bgpix) ; bgpix+(iimg2-1)*(npix2+bgpix)+1 iimg2*(npix2+bgpix)]
            if ~isempty(imgs{iimg1,iimg2})
                img( bgpix+(iimg1-1)*(npix1+bgpix)+1 : iimg1*(npix1+bgpix) , bgpix+(iimg2-1)*(npix2+bgpix)+1 : iimg2*(npix2+bgpix) ) = imgs{iimg1,iimg2};
                alp( bgpix+(iimg1-1)*(npix1+bgpix)+1 : iimg1*(npix1+bgpix) , bgpix+(iimg2-1)*(npix2+bgpix)+1 : iimg2*(npix2+bgpix) ) = alps{iimg1,iimg2};
            else
                img( bgpix+(iimg1-1)*(npix1+bgpix)+1 : iimg1*(npix1+bgpix) , bgpix+(iimg2-1)*(npix2+bgpix)+1 : iimg2*(npix2+bgpix) ) = ones(npix1,npix2)*min(min(imgs{1,1}));
                alp( bgpix+(iimg1-1)*(npix1+bgpix)+1 : iimg1*(npix1+bgpix) , bgpix+(iimg2-1)*(npix2+bgpix)+1 : iimg2*(npix2+bgpix) ) = zeros(npix1,npix2);
            end
		end
	end

	if nargin < 5
		imagesc(img);
	else
		imagesc(img, 'AlphaData',alp, 'AlphaDataMapping','none');
    end

	if nargin < 2
		clim = [min(img(:)) max(img(:))];
	end
	set(gca,'CLim',clim);

	for ibgpix=1:bgpix
		for iimg1=0:nimg1
			line(1:size(img,2), (iimg1*(npix1+bgpix)+ibgpix)*ones(1,size(img,2)), 'color',bgclr);
		end
		for iimg2=0:nimg2
			line((iimg2*(npix2+bgpix)+ibgpix)*ones(1,size(img,1)), 1:size(img,1), 'color',bgclr);
		end
	end
