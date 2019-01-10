%% plot box

function ret = PlotBox(x, y, clr, width)

if nargin < 3
	clr = [0 0 0];
end
if nargin < 4
	width = 1;
end
	
	nx = length(x);
	ny = length(y);

	line(x,y(1)*ones(1,nx),'color',clr, 'LineWidth',width);
	line(x,y(end)*ones(1,nx),'color',clr, 'LineWidth',width);
	line(x(1)*ones(1,ny),y,'color',clr, 'LineWidth',width);
	line(x(end)*ones(1,ny),y,'color',clr, 'LineWidth',width);
