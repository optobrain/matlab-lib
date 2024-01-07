clear;  close all;  clc;
comInit;

%% Clustered data with normally distributed residuals
%   We can use fitlme in this case, but let's confirm BootLME produces similar results.


% Create fake data to simulate how Math score is affected by IQ and height, when measured 
%   school by school, and how the relation differs between states.

% set parameters
clear;
states = ["A", "B", "C"];
nSch = [5 6 7];  % number of recruited schools for each state
nStud = 10;  % number of students for each school
schRndSig = 3;  % SD of the school-specific random effect
residSig = 1;  % SD of the random noise e (residual) 

% true values to be revealed by analysis
bState = [0 schRndSig/2 schRndSig*1.5];  % bonus by each state, state B weak, state C stronger
bIq = .5;  % the slope by IQ
bHeight = 0;  % bonus by height (no effect)

% create a data table
rng('default');  % for reproducibility
tb = table;
for ist=1:numel(states)  % for each state
    
    bSch = normrnd(0, schRndSig, nSch(ist), 1);  % random bonus by school (school-specific random effects)
    for isch=1:nSch(ist)  % for each school
        
        % random IQ and height of students in this school: mean IQ and height do not
        % differ between schools and states
        iq1 = normrnd(100, 10, nStud, 1);  
        height1 = normrnd(150, 20, nStud, 1);
        
        % make a small table to be added to tb
        tb1 = table;  
        tb1.state = repmat(states(ist), [nStud 1]);
        tb1.school = repmat(sprintf("sch%d", isch), [nStud 1]);
        tb1.iq = iq1;
        tb1.height = height1;
        tb1.score = bIq*iq1 + bHeight*height1 + bSch(isch) + bState(ist);
        
        tb = [tb; tb1];
        
    end
end
tb.score = tb.score + normrnd(0, residSig, size(tb,1), 1);  % add normally distributed residuals
tb.state = categorical(tb.state);
tb.school = categorical(tb.school);


% See how LM fails

% fit with all predictor variables
lm = fitlm(tb, 'score ~ iq + height + state')

% remove outliers
% two methods are widely used: Cook's distance and Leverage. 
% use the one producing less numbes of outliers (for more conservative and inclusive
% analysis)
NewFig2(2,4);
subplot(121);
    lm.plotDiagnostics('cookd');
    ioOut = find(lm.Diagnostics.CooksDistance > 3*mean(lm.Diagnostics.CooksDistance));  % indices of observation of outliers
    line(ioOut, lm.Diagnostics.CooksDistance(ioOut), Marker='o', LineStyle='none')
    title(sprintf("%d outliers", numel(ioOut)));
subplot(122);
    lm.plotDiagnostics('leverage');
    ioOutL = find(lm.Diagnostics.Leverage > 2*lm.NumCoefficients/lm.NumObservations);  % indices of observation of outliers
    line(ioOutL, lm.Diagnostics.Leverage(ioOutL), Marker='o', LineStyle='none')
    title(sprintf("%d outliers", numel(ioOutL)));
if numel(ioOutL) < numel(ioOut)
    ioOut = ioOutL;
    disp("Leverage is used instead of Cook's distance.");
end
lm1 = fitlm(tb, 'score ~ iq + height + state', Exclude=ioOut)
disp([lm.Coefficients.Estimate, lm1.Coefficients.Estimate])
disp([lm.Coefficients.pValue, lm1.Coefficients.pValue])

% height seems insignificant. 
%   when you have more than one *insignificant, numerical* predictor variables, you may 
%   want to test if any of them are correlated to each other. 
%   when correlated, you may want to remove those correlated to a key one(s).
%   see an example in: https://www.mathworks.com/help/stats/mnrfit.html#btqjoj2-1

% test if we can remove height as it has p > 0.05
if lm1.coefTest([0 0 0 0 1]) < 0.05
    lm2 = lm1;
else
    lm2 = lm1.removeTerms("height")
end
% turns out height is insignificant. remove it.

% test if we can remove school as one of them has p > 0.05
if lm2.coefTest([0 1 0 0; 0 0 1 0]) < 0.05
    lm3 = lm2;
else
    lm3 = lm2.removeTerms("state")
end
% turns out state is insignificant. remove it.

% test if the residuals are normally distributed
if ~IsNormDist(lm3.Residuals.Raw)
    warning("Residuals are not normally distributed.");
end

% compare the true vs estimated 
% here we use lm2 than lm3 to show the p value of state effect
tbRes = table;
tbRes.predictor = ["State B", "State C", "IQ"]';
tbRes.true = [bState(2:end)'; bIq];
tbRes.estimate = lm2.Coefficients.Estimate(2:end);
ci = lm2.coefCI;
tbRes.ciLow = ci(2:end,1);
tbRes.ciHigh = ci(2:end,2);
tbRes.p = lm2.Coefficients.pValue(2:end);
tbRes
% tuns out the state effect is wrong. 


% See how LME succeeds

% fit with all predictor variables
lme = fitlme(tb, 'score ~ iq + height + state + (1|school)')

% remove outliers. fitlme does not provide Diagnostics property
[~, bOut] = rmoutliers(lme.residuals);
ioOut = find(bOut);  % indices of observation of outliers
NewFig2(2,2);
line(lme.fitted, lme.residuals, Marker='x', LineStyle='none');
if sum(bOut) > 0
    line(lme.fitted(bOut), lme.residauls(bOut), Marker='o', Color='r', LineStyle='none')
end
ax = gca;
ax.YLim = [-1 1]*max(abs(ax.YLim));
line(ax.XLim, [0 0], Color='k');
ax.XLabel.String = "fitted";
ax.YLabel.String = "residuals";
ax.Title.String = sprintf("%d outliers detected", sum(bOut));
if sum(bOut) > 0
    lme1 = fitlme(tb, 'score ~ iq + height + state + (1|school)', Exclude=ioOut)
else
    lme1 = lme;
end

% height seems insignificant
%   test if any insignificant numerical predictor variables are correlated with each other
%   (no need here as we have only one insignificant numerical predictor)

% test if we can remove height as it has p > 0.05
if lme1.coefTest([0 0 0 0 1]) < 0.05
    lme2 = lme1;
else
    lme2 = fitlme(tb, 'score ~ iq + state + (1|school)', Exclude=ioOut)
end
% turns out height is insignificant. remove it.

% test if we can remove state as one of them has p > 0.05
if lme2.coefTest([0 1 0 0; 0 0 1 0]) < 0.05
    lme3 = lme2;
else
    lme3 = fitlme(tb, 'score ~ iq + (1|school)', Exclude=ioOut)
end
% turns out state is significant. keep it.

% test if adding a school-specific random effect on the slope of iq improves the model
lme4 = fitlme(tb, 'score ~ iq + state + (1|school) + (iq-1|school)', Exclude=ioOut);
cmp = lme3.compare(lme4, CheckNesting=true);  % lme3 is better than lem2?
if cmp.pValue(2) < 0.05  % yes, lme3 is better than lme2
    warning("add (iq-1|school)");
else
    disp("(iq-1|school) is insignificant random effect.");
    lme4 = lme3;
end
% turns out (iq-1|school) is insignificant. remove it.

% test if adding quadratic fixed-effect terms improves the model
lme5 = fitlme(tb, 'score ~ iq + iq^2 + state + (1|school)', Exclude=ioOut);
cmp = lme4.compare(lme5, CheckNesting=true);
if cmp.pValue(2) < 0.05
    warning("add iq^2");
else
    disp("iq^2 is insignificant fixed effect.");
    lme5 = lme4;
end
% turns out iq^2 is insignificant. remove it.

% compare the true vs estimated 
tbRes = table;
tbRes.predictor = ["State B", "State C", "IQ"]';
tbRes.true = [bState(2:end)'; bIq];
tbRes.estimate = lme5.Coefficients.Estimate(2:end);
ci = lme5.coefCI;
tbRes.ciLow = ci(2:end,1);
tbRes.ciHigh = ci(2:end,2);
tbRes.p = lme5.Coefficients.pValue(2:end);
tbRes
% now reveals effect of both states are significant, although the CI does not include the
% true value in state C.


% See how BootLME works







%%




bTrue = [1 2]';  
nGrp = 2;  
nSubj = 5;  % per group
nt = 10;  % time point per subject
eSig = .5;  % SD of the random noise e (residual)

rng('default');  % for reproducibility
t = linspace(0,1,nt)';
% t = rand(nt,1);
b01 = normrnd(0, 1, 1, nSubj);  % subject-specific intercepts
b02 = normrnd(0, 1, 1, nSubj);  % subject-specific intercepts

Y = zeros(nt, nSubj, nGrp);  X = Y;  S = strings(size(Y));  G = S;
for ig=1:nGrp
    Y(:,:,ig) = bTrue(1) + repmat(b01, [nt 1]) + ig*bTrue(2) * repmat(t, [1 nSubj]) + normrnd(0, eSig, nt, nSubj);  % [nt nSubj]
    X(:,:,ig) = repmat(t, [1 nSubj]);
    S(:,:,ig) = string(repmat((1:nSubj), [nt 1]));
    G(:,:,ig) = string(ig*ones(nt,nSubj));
end

% plot
clr = lines;
fig = NewFig2(1,nGrp);
for ig=1:nGrp
    subplot(1,nGrp,ig);
    for is=1:nSubj
        line(X(:,is,ig), Y(:,is,ig), color=clr(is,:), marker='.');
    end
    ax = gca;
    ax.XLabel.String = "x";
    ax.YLabel.String = "y";
    ax.Title.String = sprintf("group %d", ig);
end

% LM, LME
% b value are identical, but p values were lower in LME

tb = table;
tb.y = Y(:);  tb.x = X(:);  tb.s = S(:);  tb.g = G(:);

lm = fitlm(tb, "y ~ 1 + g*x")
lme = fitlme(tb, "y ~ 1 + g*x + (1|s)")
% lm = fitlm(tb, "y ~ 1 + x + g:x")  % with this more informed model, LM works ok
% lme = fitlme(tb, "y ~ 1 + x + g:x + (1|s)")

% LM, boot-LM for each group
% b values are similar. p values are lower in boot-LME

nBoot = 5000;
for ig=1:nGrp
    fprintf("\n ## [bLM, pLM, bBoot, pBoot] from LM and boot-LM for group %d: \n", ig);
    Y1 = Y(:,:,ig);  X1 = X(:,:,ig);
    [blm, lm] = BootLM(nBoot, Y1(:), [ones(numel(X1),1), X1(:)]);
    [lm.Coefficients.Estimate, lm.Coefficients.pValue, blm.b, blm.p]
end

%%

nBoot = 5000;
md = "y ~ 1 + g*x + (1|s)";
varNames = {"y", ["x"], "s", "g"};  % 4-element cell including y column name in the table, x name(s), subject name, and group name
oResamp = "r";  
bResampResid = true;  % whether to resample residuals (true) or resample X and y
ciType = "norm";
ciAlpha = 0.05;
bFig = true;
[blm, lme, fig] = BootLME(nBoot, tb, md, varNames, ...  % required inputs
    oResamp='s', ciType='norm', ciAlpha=0.05, bFig=true, ...  % options for bootstrap
    CovariancePattern='Diagonal', FitMethod='REML');  % options to be passed to fitlme()
% [blm, lme, fig] = BootLME(nBoot, tb, md, varNames, oResamp, bResampResid, ciType, ciAlpha, bFig, ...
%     CovariancePattern='Diagonal', FitMethod='REML');
% [blm, lme, fig] = BootLME(nBoot, tb, md, varNames, oResamp, bResampResid, ciType, ciAlpha, bFig, 'CovariancePattern', 'Diagonal');

% blme = BootLME(nBoot, Y(:), [ones(numel(X),1), X(:)], S(:), G(:));

%% LME

tb = table;
tb.y = Y(:);  tb.x = X(:);  tb.s = S(:);
lme = fitlme(tb, "y ~ x + (1|s)");

% cluster block bootstrap




%%

lm

% Residual bootstrap LM

bResampResid = true;  % resample only residuals? true when X is fixed.
bStudResid = false;
ciType = "norm";
[blmResid, ~, ~, fig] = BootLM(nBoot, y, X, bResampResid, bStudResid, ciType, ciAlpha, bFig);

bResampResid = true;  % resample only residuals? true when X is fixed.
bStudResid = true;
ciType = "norm";
[blmResidS, ~, ~, fig] = BootLM(nBoot, y, X, bResampResid, bStudResid, ciType, ciAlpha, bFig);

% Results

fprintf("\n ## b from [true, simple-LM, bootstrap-LM, resid-boot-LM, stud-resid-boot-LM] (when eSig=%.1f): \n", eSig);
[bTrue, lm.Coefficients.Estimate, blm.b, blmResid.b, blmResidS.b]
fprintf("\n ## b error of [simple-LM, bootstrap-LM, resid-boot-LM, stud-resid-boot-LM]: \n");
[lm.Coefficients.Estimate, blm.b, blmResid.b, blmResidS.b] - repmat(bTrue, [1 4])
fprintf("\n ## p from [simple-LM, bootstrap-LM, resid-boot-LM, stud-resid-boot-LM]: \n");
[lm.Coefficients.pValue, blm.p, blmResid.p, blmResidS.p]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from simple-LM: \n");
[ciLM, prod([ciLM(:,1)<=bTrue, ciLM(:,2)>=bTrue],2), diff(ciLM,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from bootstrap-LM: \n");
[blm.bCI, prod([blm.bCI(:,1)<=bTrue, blm.bCI(:,2)>=bTrue],2), diff(blm.bCI,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from resid-boot-LM: \n");
[blmResid.bCI, prod([blmResid.bCI(:,1)<=bTrue, blmResid.bCI(:,2)>=bTrue],2), diff(blmResid.bCI,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from stud-resid-boot-LM: \n");
[blmResidS.bCI, prod([blmResidS.bCI(:,1)<=bTrue, blmResidS.bCI(:,2)>=bTrue],2), diff(blmResidS.bCI,1,2)]







%% FUNCTIONS

function [lm, blm, blmResid, blmResidS] = TestBootLM(bTrue, nObs, eSig)
    rng('shuffle');  % to produce varying results
    x1 = rand(nObs,1);
    x2 = rand(nObs,1);
%     x1 = linspace(0, 1, nObs)';  x2 = x1;  % as fixed X : caused machine precision error
    X = [ones(nObs,1), x1, x2];  % x1 and x2 are not correlated
    y = X * bTrue + normrnd(0, eSig, nObs, 1);

    % Bootstrap LM vs simple LM

    nBoot = 5000;
    bResampResid = false;  % resample only residuals? true when X is fixed.
    bStudResid = false;
    ciType = "bca";
    ciAlpha = 0.05;
    bFig = false;  % plot a figure
    [blm, lm, ciLM] = BootLM(nBoot, y, X, bResampResid, bStudResid, ciType, ciAlpha, bFig);

    bResampResid = true;  % resample only residuals? true when X is fixed.
    bStudResid = false;
    ciType = "norm";
    blmResid = BootLM(nBoot, y, X, bResampResid, bStudResid, ciType, ciAlpha, bFig);

    bResampResid = true;  % resample only residuals? true when X is fixed.
    bStudResid = true;
    ciType = "norm";
    blmResidS = BootLM(nBoot, y, X, bResampResid, bStudResid, ciType, ciAlpha, bFig);
end


function [glm, blm, blmResid, blmResidS] = TestBootLM2(bTrue, nObs, eMean)
    rng('shuffle');  % to produce varying results
    x1 = rand(nObs,1);
    x2 = rand(nObs,1);
    X = [ones(nObs,1), x1, x2];  % x1 and x2 are not correlated
    y = X * bTrue + exprnd(eMean, nObs, 1);

    % GLM as a groundtruth?

    distr = "poisson";
    linkFun = "identity";
    glm = fitglm([x1, x2], y, Distribution=distr, Link=linkFun);

    % 4 LMs: no-boot, bootstrap, residual-boot, studentized-residual-boot

    nBoot = 5000;
    bResampResid = false;  % resample only residuals? true when X is fixed.
    bStudResid = false;
    blm = BootLM(nBoot, y, X, bResampResid, bStudResid, "bca");
    blmResid = BootLM(nBoot, y, X, true, false, "norm");
    blmResidS = BootLM(nBoot, y, X, true, true, "norm");
end
