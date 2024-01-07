%% Data with non-normally distributed residuals

clear;  close all;  clc;
comInit;


%% Create fake data to simulate how Math score is affected by IQ, height, weight, and money, 
%   and how this relationship differs between states.
%   Here the score is affected by height, which is just for simulation purpose, not true

% set parameters
states = ["A", "B", "C"];
nStud  = [10 15 20];  % number of recruited students for each state
residMean = 1;  % Mean of the exponential random noise e (residual) 

% true values to be revealed by analysis
bState = [0 residMean*0.5 residMean*1.5];  % bonus by each state, state B weak, state C stronger
bIq = .5;  % the slope of IQ
bHeight = .2;  % the slope of height

% two correlated predictor variables (we must identify them and remove one)
cWeightHeight = 0.5;  % the slope between correlated weight and height: height = cWeightHeight * height + e

% create a data table
rng('default');  % for reproducibility
% rng('shuffle');  % when want to test for different noise
tb = table;
for ist=1:numel(states)  % for each state
    
    % random IQ and height of students in this state: mean IQ and height do not
    % differ between states
    iq1 = normrnd(100, 10, nStud(ist), 1);  
    height1 = normrnd(150, 20, nStud(ist), 1);
    weight1 = cWeightHeight * height1 + normrnd(0, 10, nStud(ist), 1);
    money1 = normrnd(100, 20, nStud(ist), 1);

    % make a small table to be added to tb
    tb1 = table;  
    tb1.state = repmat(states(ist), [nStud(ist) 1]);
    tb1.iq = iq1;
    tb1.height = height1;
    tb1.weight = weight1;
    tb1.money = money1;
    tb1.score = bIq*iq1 + bHeight*height1 + bState(ist);

    tb = [tb; tb1];
end
% tb.score = tb.score + exprnd(residMean, size(tb,1), 1);  % add normally distributed residuals
tb.score = poissrnd(tb.score-min(tb.score));
tb.state = categorical(tb.state);
head(tb)


%% See how LM fails

% fit with all predictor variables
lm = fitlm(tb, 'score ~ iq + height + weight + money + state')

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
if isempty(ioOut)
    lm1 = lm;
else
    lm1 = fitlm(tb, 'score ~ iq + height + weight + money + state', Exclude=ioOut)
    disp([lm.Coefficients.Estimate, lm1.Coefficients.Estimate])
    disp([lm.Coefficients.pValue, lm1.Coefficients.pValue])
end

% weight and money seems insignificant (p>0.05). test correlation
%   when you have more than one insignificant predictor variables, you may 
%   want to test if any of them are correlated to each other. 
%   when two predictor variables are correlated, you may want to keep only one of them
%   I anticipated both height and weight are insignificant, so we can reveal their
%   correlation, but fitlm already finds height is significant.
%   see an example in: https://www.mathworks.com/help/stats/mnrfit.html#btqjoj2-1
NewFig2(2,2);
scatter(tb.weight, tb.money);
ax = gca;
ax.XLabel.String = "Weight";
ax.YLabel.String = "Money";
[r, p] = corrcoef(tb.weight, tb.money);
if p(1,2) < 0.05
    warning("Weight and money are correlated. Consider remove one of them when possible");
else
    disp("Weight and money are not correlated.");
end

% as no pair between insigificant predictors are correlated, now test if we can remove
% those insignificant one by one
% test if we can remove weight as it has p > 0.05
if lm1.coefTest([0 0 0 0 0 1 0]) < 0.05
    lm2 = lm1;
else
    lm2 = lm1.removeTerms("weight")
end
% turns out weight is insignificant, so removed it.

% test if we can remove money as it has p > 0.05
if lm2.coefTest([0 0 0 0 0 1]) < 0.05
    lm3 = lm2;
else
    lm3 = lm2.removeTerms("money")
end
% turns out money is insignificant, so removed it.

% check if adding a quadratic terms improves the model
lm4 = fitlm(tb, "score ~ state + iq + height + iq^2 + height^2", Exclude=ioOut);
if lm4.coefTest([0 0 0 0 0 1 0]) < 0.05 && lm4.coefTest([0 0 0 0 0 0 1]) < 0.05
    disp("Both quadratic terms are significant. Keep them.");
    lm5 = lm4;
elseif lm4.coefTest([0 0 0 0 0 1 0]) < 0.05 
    disp("Only iq^2 is significant. Remove height^2");
    lm5 = fitlm(tb, "score ~ state + iq + height + iq^2", Exclude=ioOut);
elseif lm4.coefTest([0 0 0 0 0 0 1]) < 0.05 
    disp("Only height^2 is significant. Remove height^2");
    lm5 = fitlm(tb, "score ~ state + iq + height + height^2", Exclude=ioOut);
else
    disp("Both quadratic terms are insignificant. Choose the previous one without them.");
    lm5 = lm3;
end

% check if the residuals are normally distributed
if ~IsNormDist(lm5.Residuals.Raw)
    warning("Residuals are not normally distributed.");
end
% turns out no (bad)

% compare the true vs estimated 
tbRes = table;
tbRes.predictor = lm5.CoefficientNames(2:end)';
trueVal = [bState(2:end)'; bIq; bHeight; 0; 0];  % including iq^2 and height^2
tbRes.true = trueVal(1:size(tbRes,1));
tbRes.estimate = lm5.Coefficients.Estimate(2:end);
ci = lm5.coefCI;
tbRes.ciLow = ci(2:end,1);
tbRes.ciHigh = ci(2:end,2);
tbRes.p = lm5.Coefficients.pValue(2:end);
tbRes.correct = ((tbRes.ciLow-tbRes.true).*(tbRes.ciHigh-tbRes.true) < 0) & (tbRes.p < 0.05);
tbResLM = tbRes
% compared to the example with normally distributed residuals (egBootLM_1.m), the
% poisson distribution made LM correctly estimate only height (wrong for all others).


%% See how GLM works when we know the true distribution of data

% fit with all predictor variables
glm = fitglm(tb, 'score ~ iq + height + weight + money + state', Distribution='poisson')

% remove outliers
% two methods are widely used: Cook's distance and Leverage. 
% use the one producing less numbes of outliers (for more conservative and inclusive
% analysis)
NewFig2(2,4);
subplot(121);
    glm.plotDiagnostics('cookd');
    ioOut = find(glm.Diagnostics.CooksDistance > 3*mean(glm.Diagnostics.CooksDistance));  % indices of observation of outliers
    line(ioOut, glm.Diagnostics.CooksDistance(ioOut), Marker='o', LineStyle='none')
    title(sprintf("%d outliers", numel(ioOut)));
subplot(122);
    glm.plotDiagnostics('leverage');
    ioOutL = find(glm.Diagnostics.Leverage > 2*glm.NumCoefficients/glm.NumObservations);  % indices of observation of outliers
    line(ioOutL, glm.Diagnostics.Leverage(ioOutL), Marker='o', LineStyle='none')
    title(sprintf("%d outliers", numel(ioOutL)));
if numel(ioOutL) < numel(ioOut)
    ioOut = ioOutL;
    disp("Leverage is used instead of Cook's distance.");
end
if isempty(ioOut)
    glm1 = glm;
else
    glm1 = fitglm(tb, 'score ~ iq + height + weight + money + state', Distribution='poisson', Exclude=ioOut)
    disp([glm.Coefficients.Estimate, glm1.Coefficients.Estimate])
    disp([glm.Coefficients.pValue, glm1.Coefficients.pValue])
end

% weight and money seems insignificant (p>0.05). test correlation
%   when you have more than one insignificant predictor variables, you may 
%   want to test if any of them are correlated to each other. 
%   when two predictor variables are correlated, you may want to keep only one of them
%   I anticipated both height and weight are insignificant, so we can reveal their
%   correlation, but fitlm already finds height is significant.
%   see an example in: https://www.mathworks.com/help/stats/mnrfit.html#btqjoj2-1
NewFig2(2,2);
scatter(tb.weight, tb.money);
ax = gca;
ax.XLabel.String = "Weight";
ax.YLabel.String = "Money";
[r, p] = corrcoef(tb.weight, tb.money);
if p(1,2) < 0.05
    warning("Weight and money are correlated. Consider remove one of them when possible");
else
    disp("Weight and money are not correlated.");
end

% as no pair between insigificant predictors are correlated, now test if we can remove
% those insignificant one by one
% test if we can remove weight as it has p > 0.05
if glm1.coefTest([0 0 0 0 0 1 0]) < 0.05
    glm2 = glm1;
else
    glm2 = glm1.removeTerms("weight")
end
% turns out weight is insignificant, so removed it.

% test if we can remove money as it has p > 0.05
if glm2.coefTest([0 0 0 0 0 1]) < 0.05
    glm3 = glm2;
else
    glm3 = glm2.removeTerms("money")
end
% turns out money is insignificant, so removed it.

% check if adding a quadratic terms improves the model
glm4 = fitglm(tb, "score ~ state + iq + height + iq^2 + height^2", Distribution='poisson', Exclude=ioOut);
if glm4.coefTest([0 0 0 0 0 1 0]) < 0.05 && glm4.coefTest([0 0 0 0 0 0 1]) < 0.05
    disp("Both quadratic terms are significant. Keep them.");
    glm5 = glm4;
elseif glm4.coefTest([0 0 0 0 0 1 0]) < 0.05 
    disp("Only iq^2 is significant. Remove height^2");
    glm5 = fitglm(tb, "score ~ state + iq + height + iq^2", Distribution='poisson', Exclude=ioOut);
elseif glm4.coefTest([0 0 0 0 0 0 1]) < 0.05 
    disp("Only height^2 is significant. Remove height^2");
    glm5 = fitglm(tb, "score ~ state + iq + height + height^2", Distribution='poisson', Exclude=ioOut);
else
    disp("Both quadratic terms are insignificant. Choose the previou sone without them.");
    glm5 = glm3;
end
% iq^2 is significant (wrong) although its estimate is very small (-0.001)

% with GLM, it's normal for residuals to non-normally distributed, so no check

% compare the true vs estimated 
tbRes = table;
tbRes.predictor = glm5.CoefficientNames(2:end)';
trueVal = [bState(2:end)'; bIq; bHeight; 0; 0];  % including iq^2 and height^2
tbRes.true = trueVal(1:size(tbRes,1));
tbRes.estimate = glm5.Coefficients.Estimate(2:end);
ci = glm5.coefCI;
tbRes.ciLow = ci(2:end,1);
tbRes.ciHigh = ci(2:end,2);
tbRes.p = glm5.Coefficients.pValue(2:end);
tbRes.correct = ((tbRes.ciLow-tbRes.true).*(tbRes.ciHigh-tbRes.true) < 0) & (tbRes.p < 0.05);
tbResGLM = tbRes
% compared to the LM result, GLM with the known distribution revealed state C effect,
% but made all estimations wrong (CIs do not include the true values).
% Although it works better in some aspects, WE CAN USE GLM ONLY WHEN WE KNOW THE DISTRIBUTION. 
% In this example, the use of GLM with the poisson distribution can be justified
% only when the math score can be considered as a "count" (e.g., how many problems a student
% correctly answer to).


%% Check if BootLM works better

% === input ===
% tb : data table to be inputted to fitlm()
% md : formula to be inputted to fitlm()
% varNames : 2-element string array:
%               varNames(1) : y variable name
%               varNames(2) : group variable name. When groups are determined by
%                             more than one variable, make a column of a unique
%                             group identifier and indicate it here. 
%               e.g,. varNames = ["y", "grpId"]

% === option === 
% oResamp : How to resample.    
%       'resid'     = [default] Resample all residuals, randomly across groups.
%                     Good when X is fixed (every predictor variable has a designed set of *fixed* values)
%       'residGrp'  = Resample residuals within each group.
%                     Not recommended unless residuals are extremely different between groups.
%       'obsGrp'    = Resample observations (original y and x values) within each group. 
%                     Good when X is random (subjects are randomly selected within group)
% bStudResid : whether to use studentized residuals. See https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Resampling_residuals
%               Default=true following the wikipedia.
% ciType : only "norm" or "per"
% ciAlpha : 0.05 (default)
% weight : The option Weights of bootstrp
% bPara : whether to use parallel computing
% bFig : whether to make a figure (distribution of bootstrapped results)

% === output ====
% blm : struct of results
%   .b : average of bootstrapped coefficients (see lme for coefficient order)
%   .ci : CI of bootstrapped coefficients
%   .p : p value of bootstrapped coefficients against null hypothesis-bootstrapped ones.
%        Not provided when oResamp='obsGrp' as it resamples *residuals*.
%   .bBoot : bootstrapped coefficient values [NumCoefficients nBoot]
%   .bBootNull : null hypothesis-bootstrapped coefficient values [NumCoefficients nBoot]
%                Not provided when oResamp='obsGrp' as it resamples *residuals*.
% lm : LM object, return of simple fitlm
% fig : figure handle

% use the best model found above with fitlme
nBoot = 2000;  % for p value resolution of 0.001
md = "score ~ state + iq + height + weight + money";
varNames = ["score", "state"];  % y, group

% bootstrap coeffiicents by resampling residuals, with default options
[blm, lm, fig] = BootLM(nBoot, tb, md, varNames, ...  % required inputs
    bPara=true, bFig=true, ...  % options for bootstrap
    Exclude=ioOut);  % options for fitlm
tbRes = table;
tbRes.predictor = ["State B", "State C", "IQ", "Height", "Weight", "Money"]';
tbRes.true = [bState(2:end)'; bIq; bHeight; 0; 0];
tbRes.estimate = blm.b(2:end);  tbRes.ciLow = blm.ci(2:end,1);  tbRes.ciHigh = blm.ci(2:end,2);  tbRes.p = blm.p(2:end);
tbRes.correct = ((tbRes.ciLow-tbRes.true).*(tbRes.ciHigh-tbRes.true) < 0) & (tbRes.p < 0.05);
tbResBoot = tbRes
% not work well, because the default oResamp='resid' is good when X is fixed

% bootstrap coeffiicents by resampling residuals within each group
[blm1, ~, fig1] = BootLM(nBoot, tb, md, varNames, ...  % required inputs
    oResamp='residGrp', bPara=true, bFig=true, ...  % options for bootstrap
    Exclude=ioOut);  % options for fitlm
tbRes.estimate = blm1.b(2:end);  tbRes.ciLow = blm1.ci(2:end,1);  tbRes.ciHigh = blm1.ci(2:end,2);  tbRes.p = blm1.p(2:end);
tbRes.correct = ((tbRes.ciLow-tbRes.true).*(tbRes.ciHigh-tbRes.true) < 0) & (tbRes.p < 0.05);
tbResBoot_residGrp = tbRes
% not work well. oResamp='residGrp' is not recommended

% bootstrap coeffiicents by resampling observations within each group
[blm2, ~, fig2] = BootLM(nBoot, tb, md, varNames, ...  % required inputs
    oResamp='obsGrp', bPara=true, bFig=true, ...  % options for bootstrap
    Exclude=ioOut);  % options for fitlm
tbRes.estimate = blm2.b(2:end);  tbRes.ciLow = blm2.ci(2:end,1);  tbRes.ciHigh = blm2.ci(2:end,2);  tbRes.p(:) = NaN;
tbRes.correct = ((tbRes.ciLow-tbRes.true).*(tbRes.ciHigh-tbRes.true) < 0) & ...
    ((abs(tbRes.true)>0 & tbRes.ciHigh.*tbRes.ciLow>0) | (tbRes.true==0 & tbRes.ciHigh.*tbRes.ciLow<0));
tbResBoot_obsGrp = tbRes
% works well

% compare all results
tbResLM
tbResGLM
tbResBoot
tbResBoot_residGrp
tbResBoot_obsGrp
% when testing with rng(shuffle), this approach produced the best result in most cases.
% In most cases, this made correct estimation for all variables except
% for state B (its effect was weaker than the noise level in the original data).


%% ----- OLD CODE
%{
clear;  close all;  clc;
comInit;


%% Simulated data with normally distributed residuals

% Create and plot data for y = b0 + b1*x1 + b2*x2 + e

bTrue = [1 2 3]';  % b0 b1 b2
nObs = 20;
eSig = 1;  % SD of the random noise e (residual)

rng('default');  % for reproducibility
x1 = rand(nObs,1);
x2 = rand(nObs,1);
% x1 = linspace(.1, .9, nObs)';  x2 = x1 + 0.05;  % as fixed X : caused machine precision error
X = [ones(nObs,1), x1, x2];  % x1 and x2 are not correlated
y = X * bTrue + normrnd(0, eSig, nObs, 1);

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

% Simple LM, bootstrap LM

nBoot = 5000;
bResampResid = false;  % resample only residuals? true when X is fixed.
bStudResid = false;
ciType = "bca";
ciAlpha = 0.05;
bFig = true;  % plot a figure
[blm, lm, ciLM, fig] = BootLM(nBoot, y, X, bResampResid, bStudResid, ciType, ciAlpha, bFig);

lm
if prod(TestNormality(lm.Residuals.Raw)),  warning("residuals are not normally distributed.");  end

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


%% compare simple LM vs bootstrap LM vs residual bootstrap LM

nTest = 100;
bErr = zeros(nTest,4);  % RMS of 3 b's
nGoodCI = zeros(nTest,4);  % how many b CIs (total 3) include the true value?
ciWid = zeros(nTest,4);  % CI width (mean of 3 CIs)

% Repeat

DispProg(0, nTest, "Repeating ...");
for it=1:nTest
    [lm, blm, blmResid, blmResidS] = TestBootLM(bTrue, nObs, eSig);
    ciLM = coefCI(lm);
    bErr(it,:) = sqrt(mean(([lm.Coefficients.Estimate, blm.b, blmResid.b, blmResidS.b] - repmat(bTrue, [1 4])).^2,1));
    nGoodCI(it,:) = [sum(prod([ciLM(:,1)<=bTrue, ciLM(:,2)>=bTrue],2)), ...
        sum(prod([blm.bCI(:,1)<=bTrue, blm.bCI(:,2)>=bTrue],2)), ...
        sum(prod([blmResid.bCI(:,1)<=bTrue, blmResid.bCI(:,2)>=bTrue],2)), ...
        sum(prod([blmResidS.bCI(:,1)<=bTrue, blmResidS.bCI(:,2)>=bTrue],2)), ...
        ];
    ciWid(it,:) = mean([diff(ciLM,1,2), diff(blm.bCI,1,2), diff(blmResid.bCI,1,2), diff(blmResidS.bCI,1,2)],1);
    if mod(it,10) == 0
        DispProg(it, nTest);
    end
end

% Results

mtit = ["LM", "boot LM", "resid-boot LM", "stud-resid-boot LM"];
xl = strings(3,4);  % [imetric im]
fig = NewFig2(2,5);
for im=1:4
    subplot(131);  histogram(bErr(:,im));  hold on;  xl(1,im) = sprintf("%.2g +/- %.2g", mean(bErr(:,im)), std(bErr(:,im)));
    subplot(132);  histogram(nGoodCI(:,im));  hold on;  xl(2,im) = sprintf("%.2g +/- %.2g", mean(nGoodCI(:,im)), std(nGoodCI(:,im)));
    subplot(133);  histogram(ciWid(:,im));  hold on;  xl(3,im) = sprintf("%.2g +/- %.2g", mean(ciWid(:,im)), std(ciWid(:,im)));
end
subplot(131);  ax = gca;  ax.XLabel.String = xl(1,:);  ax.Title.String = "b error RMS";
subplot(132);  ax = gca;  ax.XLabel.String = xl(2,:);  ax.Title.String = "Num good CI";  legend(mtit, location="northwest");
subplot(133);  ax = gca;  ax.XLabel.String = xl(3,:);  ax.Title.String = "CI width";

fprintf("\n ## boot-LM vs resid-boot-LM: \n");
fprintf("b error RMS: p=%.2g \n", BootCompareTwo(bErr(:,2), bErr(:,3), nBoot));
fprintf("Num good CI: p=%.2g \n", BootCompareTwo(nGoodCI(:,2), nGoodCI(:,3), nBoot));
fprintf("CI width: p=%.2g \n", BootCompareTwo(ciWid(:,2), ciWid(:,3), nBoot));

fprintf("\n ## resid-boot-LM vs stud-resid-boot-LM: \n");
fprintf("b error RMS: p=%.2g \n", BootCompareTwo(bErr(:,3), bErr(:,4), nBoot));
fprintf("Num good CI: p=%.2g \n", BootCompareTwo(nGoodCI(:,3), nGoodCI(:,4), nBoot));
fprintf("CI width: p=%.2g \n", BootCompareTwo(ciWid(:,3), ciWid(:,4), nBoot));

% Findings
%   The b error was not different across three bootstrap LMs. Between resid-boot vs
%   stud-resid-boot, the stud-resid-boot produced significantly narrower CIs but slightly
%   lower good CI numbers. Stud-resid-boot looks better.
%   Between boot vs stud-resid-boot, you must choose stud-resid-boot-LM when X is fixed 
%   (controlled experiment, like a set of evenly-distributed drug dose values), 
%   while choosing boot-LM when X is random (sampling from population, like height [x] vs weight [y] of people).


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

% GLM as a groundtruth?

distr = "poisson";
linkFun = "identity";
glm = fitglm([x1, x2], y, Distribution=distr, Link=linkFun)
ci = glm.coefCI;
if prod(TestNormality(glm.Residuals.Raw)),  warning("residuals of glm are not normally distributed.");  end

% 4 LMs: no-boot, bootstrap, residual-boot, studentized-residual-boot

nBoot = 5000;
bResampResid = false;  % resample only residuals? true when X is fixed.
bStudResid = false;
[blm, lm, ciLM] = BootLM(nBoot, y, X, bResampResid, bStudResid, "bca");
blmResid = BootLM(nBoot, y, X, true, false, "norm");
blmResidS = BootLM(nBoot, y, X, true, true, "norm");
if prod(TestNormality(lm.Residuals.Raw)),  warning("residuals of lm are not normally distributed.");  end

% Results

fprintf("\n ## b from [true, GLM, simple-LM, bootstrap-LM, resid-boot-LM, stud-resid-boot-LM] (when eMean=%.1f): \n", eMean);
[bTrue, glm.Coefficients.Estimate, lm.Coefficients.Estimate, blm.b, blmResid.b, blmResidS.b]
fprintf("\n ## b error of [GLM, simple-LM, bootstrap-LM, resid-boot-LM, stud-resid-boot-LM]: \n");
[glm.Coefficients.Estimate, lm.Coefficients.Estimate, blm.b, blmResid.b, blmResidS.b] - repmat(bTrue, [1 5])
fprintf("\n ## p from [GLM, simple-LM, bootstrap-LM, resid-boot-LM, stud-resid-boot-LM]: \n");
[glm.Coefficients.pValue, lm.Coefficients.pValue, blm.p, blmResid.p, blmResidS.p]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from GLM: \n");
[ci, prod([ci(:,1)<=bTrue, ci(:,2)>=bTrue],2), diff(ci,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from simple-LM: \n");
[ciLM, prod([ciLM(:,1)<=bTrue, ciLM(:,2)>=bTrue],2), diff(ciLM,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from bootstrap-LM: \n");
[blm.bCI, prod([blm.bCI(:,1)<=bTrue, blm.bCI(:,2)>=bTrue],2), diff(blm.bCI,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from resid-boot-LM: \n");
[blmResid.bCI, prod([blmResid.bCI(:,1)<=bTrue, blmResid.bCI(:,2)>=bTrue],2), diff(blmResid.bCI,1,2)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from stud-resid-boot-LM: \n");
[blmResidS.bCI, prod([blmResidS.bCI(:,1)<=bTrue, blmResidS.bCI(:,2)>=bTrue],2), diff(blmResidS.bCI,1,2)]


%% compare GLM vs bootstrap LM vs residual bootstrap LM

nTest = 100;
bErr = zeros(nTest,4);  % RMS of 3 b's
nGoodCI = zeros(nTest,4);  % how many b CIs (total 3) include the true value?
ciWid = zeros(nTest,4);  % CI width (mean of 3 CIs)

% Repeat

DispProg(0, nTest, "Repeating ...");
for it=1:nTest
    [glm, blm, blmResid, blmResidS] = TestBootLM2(bTrue, nObs, eMean);
    ci = glm.coefCI;
    bErr(it,:) = sqrt(mean(([glm.Coefficients.Estimate, blm.b, blmResid.b, blmResidS.b] - repmat(bTrue, [1 4])).^2,1));
    nGoodCI(it,:) = [sum(prod([ci(:,1)<=bTrue, ci(:,2)>=bTrue],2)), ...
        sum(prod([blm.bCI(:,1)<=bTrue, blm.bCI(:,2)>=bTrue],2)), ...
        sum(prod([blmResid.bCI(:,1)<=bTrue, blmResid.bCI(:,2)>=bTrue],2)), ...
        sum(prod([blmResidS.bCI(:,1)<=bTrue, blmResidS.bCI(:,2)>=bTrue],2)), ...
        ];
    ciWid(it,:) = mean([diff(ci,1,2), diff(blm.bCI,1,2), diff(blmResid.bCI,1,2), diff(blmResidS.bCI,1,2)],1);
    if mod(it,10) == 0
        DispProg(it, nTest);
    end
end

% Results

mtit = ["GLM", "boot LM", "resid-boot LM", "stud-resid-boot LM"];
xl = strings(3,4);  % [imetric im]
fig = NewFig2(2,5);
for im=1:4
    subplot(131);  histogram(bErr(:,im));  hold on;  xl(1,im) = sprintf("%.2g +/- %.2g", mean(bErr(:,im)), std(bErr(:,im)));
    subplot(132);  histogram(nGoodCI(:,im));  hold on;  xl(2,im) = sprintf("%.2g +/- %.2g", mean(nGoodCI(:,im)), std(nGoodCI(:,im)));
    subplot(133);  histogram(ciWid(:,im));  hold on;  xl(3,im) = sprintf("%.2g +/- %.2g", mean(ciWid(:,im)), std(ciWid(:,im)));
end
subplot(131);  ax = gca;  ax.XLabel.String = xl(1,:);  ax.Title.String = "b error RMS";
subplot(132);  ax = gca;  ax.XLabel.String = xl(2,:);  ax.Title.String = "Num good CI";  legend(mtit, location="northwest");
subplot(133);  ax = gca;  ax.XLabel.String = xl(3,:);  ax.Title.String = "CI width";

fprintf("\n ## boot-LM vs resid-boot-LM: \n");
fprintf("b error RMS: p=%.2g \n", BootCompareTwo(bErr(:,2), bErr(:,3), nBoot));
fprintf("Num good CI: p=%.2g \n", BootCompareTwo(nGoodCI(:,2), nGoodCI(:,3), nBoot));
fprintf("CI width: p=%.2g \n", BootCompareTwo(ciWid(:,2), ciWid(:,3), nBoot));

fprintf("\n ## resid-boot-LM vs stud-resid-boot-LM: \n");
fprintf("b error RMS: p=%.2g \n", BootCompareTwo(bErr(:,3), bErr(:,4), nBoot));
fprintf("Num good CI: p=%.2g \n", BootCompareTwo(nGoodCI(:,3), nGoodCI(:,4), nBoot));
fprintf("CI width: p=%.2g \n", BootCompareTwo(ciWid(:,3), ciWid(:,4), nBoot));

% Findings
%   The b error was not different across three bootstrap LMs. Between resid-boot vs
%   stud-resid-boot, the stud-resid-boot produced significantly broader CIs (still
%   narrower than GLM and XY-boot-LM), but it produced significantly higher good CI numbers. 
%   Stud-resid-boot looks better.
%   Between boot vs stud-resid-boot, stud-resid-boot produced slightly narrower CIs, but 
%   you must choose one based on the nature of data: stud-resid-boot-LM when X is fixed 
%   (controlled experiment, like a set of evenly-distributed drug dose values), 
%   boot-LM when X is random (sampling from population, like height [x] vs weight [y] of people).

% Default setting
%   Based on the above findings, BootLM()'s default setting is bReampResid=true and
%   bStudResid=true.


%% Simulated data clustered by subjects with normally distributed residuals

clear;

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
%}
