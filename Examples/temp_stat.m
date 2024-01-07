clear;

load fisheriris

X = meas(51:end,:);
y = strcmp('versicolor', species(51:end));

glm = fitglm(X, y, 'linear', Distribution='binomial')

glm.coefCI
p = coefTest(glm)

p = coefTest(glm, [0 0 0 1 0])
p = coefTest(glm, [0 1 0 0 0])
p = coefTest(glm, [0 0 1 0 0])
p = coefTest(glm, [0 1 0 0 0; 0 0 1 0 0])
 
%%

clear

x = sort(rand(20,1));
yi = sort(rand(20,1));
catnames = ["small", "medium", "large"];
y = discretize(yi, [0 .4 .7 1], 'categorical',catnames);  % ordinal category

fitglm(x, y, Distribution='binomial')


%% 

clear

rng('default') % for reproducibility
X = randn(100,7);
mu = exp(X(:,[1 3 6])*[.4;.2;.3] + 1);  % E(y) is affected by x1, x3, and x6 with slopes of .4, .2, and .3
y = poissrnd(mu);
%  Fit a generalized linear model using the Poisson distribution.

glm = fitglm(X,y,'linear','Distribution','poisson')
lm = fitlm(X,y,'linear')