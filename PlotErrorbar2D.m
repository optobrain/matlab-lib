% nbin > 0 :: split data by nbin by the order of size
% x0 and x00 :: for idx for x0 = x00(ix)

function [Ym Xm Yd Xd] = PlotErrorbar2D(x, y, nbin, x0, x00, clr, mkr)

if nargin < 7
	mkr = '.';
end
if nargin < 6
	clr = 'b';
end
if nargin < 5
	x00 = [];
else
    x00 = x00(:);
end
if nargin < 4
	x0 = [];
end
if nargin < 3
	nbin = 0;
end

if size(x) ~= size(y)
	error('size mismatch');
end
if numel(x0) > 0
	if size(x) ~= size(x0)
		error('size mismatch');
	end
end

	x = x(:);  nx = length(x);  
	y = y(:);
	if nbin > 0
		[x is] = sort(x);  y = y(is);
		nxb = floor(nx/nbin);
		Xm = zeros(1,nbin);  Xd = Xm;  Ym = Xm;  Yd = Xm;
		for ib=1:nbin
			ix = (ib-1)*nxb+1:ib*nxb;  
			xx = x(ix);  yy = y(ix); 
			xm = mean(xx);  xd = std(xx);  ym = mean(yy);  yd = std(yy);
			line(xm,ym,'Marker',mkr,'color',clr);
			line(xm*ones(1,100),linspace(ym-yd,ym+yd,100),'color',clr);
			line(linspace(xm-xd,xm+xd,100),ym*ones(1,100),'color',clr);
			Xm(ib) = xm;  Xd(ib) = xd;  Ym(ib) = ym;  Yd(ib) = yd;
		end
	else	% use x0 and x00
		x0 = x0(:);
		if numel(x00) > 0
			Xm = zeros(1,length(x00));  Xd = Xm;  Ym = Xm;  Yd = Xm;
			for ib=1:length(x00)
				ix = find(x0 == x00(ib));
				xx = x(ix);  yy = y(ix);
				xm = mean(xx);  xd = std(xx);  ym = mean(yy);  yd = std(yy);
				line(xm,ym,'Marker',mkr,'color',clr);
				line(xm*ones(1,100),linspace(ym-yd,ym+yd,100),'color',clr);
				line(linspace(xm-xd,xm+xd,100),ym*ones(1,100),'color',clr);
				Xm(ib) = xm;  Xd(ib) = xd;  Ym(ib) = ym;  Yd(ib) = yd;
			end
		else
			% split by unique in x0
			if numel(x0) == 0
				x0 = x;
			end
			[x0 is] = sort(x0);  x = x(is);  y = y(is);
			uix = find([x0(1); diff(x0)] > 0);
			Xm = zeros(1,length(uix));  Xd = Xm;  Ym = Xm;  Yd = Xm;
			for ib=1:length(uix)
				if ib < length(uix)
					ix = uix(ib):uix(ib+1)-1;
				else
					ix = uix(ib):nx;
				end
				xx = x(ix);  yy = y(ix);
				xm = mean(xx);  xd = std(xx);  ym = mean(yy);  yd = std(yy);
				line(xm,ym,'Marker',mkr,'color',clr);
				line(xm*ones(1,100),linspace(ym-yd,ym+yd,100),'color',clr);
				line(linspace(xm-xd,xm+xd,100),ym*ones(1,100),'color',clr);
				Xm(ib) = xm;  Xd(ib) = xd;  Ym(ib) = ym;  Yd(ib) = yd;
			end
		end
	end

