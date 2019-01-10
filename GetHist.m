
% x: data
% bNorm: hy = hy / length(x)?
% nb: nbin 
% edgeFirst: the lower edge of the first bin (when nbin = [])
% edgeArray: edges of the middle bins (when nbin = [])
% edgeEnd: the higher edge of the last bin (when nbin = [])

function [hy,hx] = GetHist(x, bNorm, nb, edgeFirst, edgeArray, edgeEnd)

    if ~isempty(nb)
        [hy,edges] = histcounts(x,nb);
        hx = edges(1:end-1)+mean(diff(edges))/2;        
    else
        nb = length(edgeArray)+1;
        edges = zeros(nb+1,1);
        edges(1) = edgeFirst;
        edges(2:end-1) = edgeArray;
        edges(end) = edgeEnd;
        hy = histcounts(x,edges);
        hx = zeros(1,nb);
        hx(1:end-1) = edgeArray-mean(diff(edgeArray))/2;
        hx(end)= edgeArray(end)+mean(diff(edgeArray))/2;
    end
    
    if bNorm
        hy = hy/length(x);
    end