%% Plot waterfall of Y(ch,t)

% vline = vertical lines at t = [0 0.1], for e.g.
% nstd = multiple of std as vertical distance
% vd = vertical distance btw lines
% sclr = line color string

function ret = PlotWaterfall(Y, t, vline, nstd, vd, sclr)

if nargin < 3
    vline = [];
end
if nargin < 5 || length(vd) == 0
    if nargin < 4 || length(nstd) == 0
        nstd = 10;
    end
    vd = sqrt(mean(std(Y,1,2).^2))*nstd;
end
if nargin < 6 || length(sclr) == 0
    sclr = 'lines';
end

    [nch nt] = size(Y);
    cl = eval([sclr '(' num2str(nch) ');']);
    
%         plot(t,Y'-repmat([1:nch]*vd, [nt 1]));
        for ich=1:nch
            line(t,Y(ich,:)-ich*vd,'color',cl(ich,:));
        end
        set(gca,'YTick',[-nch:0]*vd);  set(gca,'YTickLabel',[[nch:-1:1] vd]);
        for iv=1:length(vline)
            line(vline(iv)*[1 1], get(gca,'ylim'), 'color',[1 1 1]*3/4);
        end
