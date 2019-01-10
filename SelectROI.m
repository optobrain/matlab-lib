

function [xx,yy] = SelectROI()

    h = imrect;
    roi = round(wait(h));
    delete(h);
    xx = roi(1)+(1:roi(3));
    yy = roi(2)+(1:roi(4));
