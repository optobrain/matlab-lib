%% plot trace in the complex plane

function ret = PlotTrace(y, clrmap, marker)

if nargin <3
	marker = '.';
end
if nargin < 2
	clrmap = 'jet';
end

nt = length(y);

clr = eval([ clrmap '(' num2str(nt) ')']);

for it=1:nt
	line(real(y(it)),imag(y(it)),'color',clr(it,:),'Marker',marker);
end

