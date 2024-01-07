clear;  close all;  clc;
comInit;


%% Simulated data with normally distributed residuals

% Create and plot data for a mixture of logistic and linear relations

mdfun = @(b,X) b(1)./(1+exp(-b(2)*(X(:,1)-b(3)))) + b(4)*X(:,2);
bTrue = [1 10 .5 1]';  
nObs = 40;
eSig = .1;  % SD of the random noise e (residual)

rng('default');  % for reproducibility. comment this when you want to repeat different simulations
x1 = rand(nObs,1);
x2 = rand(nObs,1);
X = [x1, x2];
y = mdfun(bTrue, X) + normrnd(0, eSig, nObs, 1);

% plot
fig = NewFig2(2,4);
subplot(121);
    plot3(x1, x2, y, linestyle='none', marker='o');
    ax = gca;
    ax.XGrid = 'on';  ax.YGrid = 'on';  ax.ZGrid = 'on';
    ax.XLabel.String = "x1";
    ax.YLabel.String = "x2";
    ax.ZLabel.String = "y";
subplot(122);
    scatter(x1, y, [], x2);
    ax = gca;
    ax.Colormap = jet;
    ax.XGrid = 'on';  ax.YGrid = 'on';  ax.ZGrid = 'on';
    ax.XLabel.String = ["x1", "color = x2"];
    ax.YLabel.String = "y";

% Bootstrap NLM vs simple NLM

nBoot = 1000;
b0 = bTrue;
ciType = "norm";
ciAlpha = 0.05;
bFig = true;  % plot a figure
[bnlm, nlm, ciNLM, fig] = BootNLM(nBoot, y, X, mdfun, b0, ciType, ciAlpha, bFig);

nlm

fprintf("\n ## b from [true, simple-LM, bootstrap-LM] (when eSig=%.1f): \n", eSig);
[bTrue, nlm.Coefficients.Estimate, bnlm.b]
fprintf("\n ## p from [simple-LM, bootstrap-LM] (when eSig=%.1f): \n", eSig);
[nlm.Coefficients.pValue, bnlm.p]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from simple-LM (when eSig=%.1f): \n", eSig);
[ciNLM, prod([ciNLM(:,1)<=bTrue, ciNLM(:,2)>=bTrue],2), diff(ciNLM,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from bootstrap-LM (when eSig=%.1f): \n", eSig);
[bnlm.bCI, prod([bnlm.bCI(:,1)<=bTrue, bnlm.bCI(:,2)>=bTrue],2), diff(bnlm.bCI,1,2)]

% findings
%   Both LM and bootstrap-LM produced good CIs that involves the true values, but the
%   bootstrap-LM produced narrower CIs, which suggests that the bootstrap-LM is better as
%   it more specifically predict the true values.


%% Simulated data with non-normally distributed residuals

% Create and plot data for y = b0 + b1*x1 + b2*x2 + e
clear;

bTrue = [1 2 3]';  % b0 b1 b2
nObs = 30;
eMean = 2;  % Mean of the exponential random noise e (residual)

rng('default');  % for reproducibility. comment this when you want to repeat different simulations
x1 = rand(nObs,1);
x2 = rand(nObs,1);
X = [ones(nObs,1), x1, x2];  % x1 and x2 are not correlated
y = X * bTrue + exprnd(eMean, nObs, 1);

% plot
fig = NewFig2(2,4);
subplot(121);
    plot3(x1, x2, y, linestyle='none', marker='o');
    ax = gca;
    ax.XGrid = 'on';  ax.YGrid = 'on';  ax.ZGrid = 'on';
    ax.XLabel.String = "x1";
    ax.YLabel.String = "x2";
    ax.ZLabel.String = "y";
subplot(122);
    scatter(x1, y, [], x2);
    ax = gca;
    ax.Colormap = jet;
    ax.XGrid = 'on';  ax.YGrid = 'on';  ax.ZGrid = 'on';
    ax.XLabel.String = ["x1", "color = x2"];
    ax.YLabel.String = "y";

% Bootstrap LM, simple LM

nBoot = 1000;
ciType = "norm";
ciAlpha = 0.05;
bFig = true;  % plot a figure
[bnlm, nlm, ciNLM, fig] = BootLM(nBoot, y, X, ciType, ciAlpha, bFig);
nlm

% Simple GLM

distr = "gamma";
linkFun = "identity";
glm = fitglm([x1, x2], y, Distribution=distr, Link=linkFun);
ciGLM = coefCI(glm);
glm

fprintf("\n ## b from [true, simple-LM, bootstrap-LM, GLM] (when eMean=%.1f): \n", eMean);
[bTrue, nlm.Coefficients.Estimate, bnlm.b, glm.Coefficients.Estimate]
fprintf("\n ## p from [simple-LM, bootstrap-LM, GLM]: \n");
[nlm.Coefficients.pValue, bnlm.p, glm.Coefficients.pValue]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from simple-LM: \n");
[ciNLM, prod([ciNLM(:,1)<=bTrue, ciNLM(:,2)>=bTrue],2), diff(ciNLM,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from bootstrap-LM: \n");
[bnlm.bCI, prod([bnlm.bCI(:,1)<=bTrue, bnlm.bCI(:,2)>=bTrue],2), diff(bnlm.bCI,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from GLM with %s-%s: \n", distr, linkFun);
[ciGLM, prod([ciGLM(:,1)<=bTrue, ciGLM(:,2)>=bTrue],2), diff(ciGLM,1,2)]

% findings
%   The bootstrap LM worked the best: its CIs involved the true values, with the narrowest
%   CIs.
%   The simple LM worked better than expected.
%   The GLM with gamma distribution worked better than the one with Poisson, but still
%   produced slightly wider CIs than the bootstrap LM. 


