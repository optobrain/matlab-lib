%% plot circle

% cor: center of rotation [x y]
% r: radius

function ret = PlotCircle(cor, r, clr, width, ntheta)

if nargin < 5
    ntheta = 20;
end
if nargin < 4
	width = 1;
end
if nargin < 3
	clr = [0 0 0];
end


    theta = linspace(0,2*pi,ntheta);
    line(cor(1)+r*cos(theta),cor(2)+r*sin(theta),'color',clr,'linewidth',width);
