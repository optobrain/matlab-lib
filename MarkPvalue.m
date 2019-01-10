
% x: x values
% y: y values (mean)
% yd: STD of y
% p: p values
% thrP: threshold p values

function ret = MarkPvalue(gca, x, y, yd, p, thrP)

if nargin < 6
    thrP = [0.05 0.01 0.001];
end
thrP = sort(thrP,'descend');

    nx = length(x);
    s = cell(nx,1);
    for ix=1:nx
        s{ix} = '';
        for it=1:length(thrP)
            if p(ix) < thrP(it)
                s{ix} = '';
                for ii=1:it
                    s{ix} = [s{ix} '*'];
                end
            end
        end
    end
    text(x,y+yd+diff(get(gca,'ylim'))/20,s,'HorizontalAlignment','center');
    ret = true;
    