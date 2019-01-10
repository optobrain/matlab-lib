% convert resnorm from lsqcurvefit into R2

function R2 = ResnormToR2(resnorm,y)

R2 = 1 - resnorm / sum((y-mean(y)).^2);