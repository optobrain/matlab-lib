clear;  close all;  clc;
comInit;


%% Simulated data of two groups with normal distribution

% Simulate data

clear;
bTrue = 0.5;
eSig = 1;
rng('default')  % for reproducibility
yCtr = normrnd(1, eSig, 200,1);  % control group
yExp = normrnd(1+bTrue, eSig, 100,1);  % experimental group

% Student t-test

fig = NewFig(1,2);
subplot(121);
    histogram(yCtr);  hold on;  histogram(yExp);
    legend(["control", "experiment"])
    [~, p] = ttest2(yCtr, yExp);
    ax = gca;
    ax.XLabel.String = [sprintf("p=%.3g", p), ...
        sprintf("E(y_{ctr})=%.2f, E(y_{exp})=%.2f", mean(yCtr), mean(yExp))];
    ax.Title.String = ["Student t-test", "on original data"];
% subplot(122);
%     histogram(log(yCtr));  hold on;  histogram(log(yExp));
%     [~, p] = ttest2(log(yCtr), log(yExp));
%     ax = gca;
%     ax.XLabel.String = sprintf("p=%.3g", p);
%     ax.Title.String = ["Student t-test", "after log transform"];

% Simple LM, Bootstrap LM with dummy var

nBoot = 2000;

y = [yCtr; yExp];
X = [zeros(size(yCtr)); ones(size(yExp))];
X = [ones(size(X)), X];
bResampResid = true;  % since X is fixed.
[blm, lm, ciLM] = BootLM(nBoot, y, X, bResampResid);  % see egBootLM.m for detailed/default options

% Bootstrap hypothesis testing

nSided = 2;  % two-sided
effectSizeType = "meanDiff";  % meanDiff, meanRatio, medianDiff, medianRatio, or linReg
effectSizeAlpha = 0.05;
bFig = true;  % plot a related figure
[pTest, t, tBoot, effectSize] = BootCompareTwo(yCtr, yExp, nBoot, nSided, effectSizeType, effectSizeAlpha, bFig);

% Results

fprintf("\n ## b (effect size) from [true, simple-LM, bootstrap-LM, bootstrap-test]: \n");
[bTrue, lm.Coefficients.Estimate(2), blm.b(2), effectSize.mean]
fprintf("\n ## b error of [simple-LM, bootstrap-LM, bootstrap-test]: \n");
[lm.Coefficients.Estimate(2), blm.b(2), effectSize.mean] - bTrue
fprintf("\n ## p from [t-test, simple-LM, bootstrap-LM, bootstrap-test]: \n");
[p, lm.Coefficients.pValue(2), blm.p(2), pTest]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from simple-LM: \n");
[ciLM(2,:), prod([ciLM(2,1)<=bTrue, ciLM(2,2)>=bTrue],2), diff(ciLM(2,:))]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from bootstrap-LM: \n");
[blm.bCI(2,:), prod([blm.bCI(2,1)<=bTrue, blm.bCI(2,2)>=bTrue],2), diff(blm.bCI(2,:))]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from bootstrap-test: \n");
ci = effectSize.ci';
[ci, prod([ci(1)<=bTrue, ci(2)>=bTrue],2), diff(ci)]


%% compare simple-LM vs bootstrap-LM vs bootstrap-test

nTest = 100;
pErr = zeros(nTest,3);  % error against t-test p in log scale
bErr = zeros(nTest,3);  
nGoodCI = zeros(nTest,3);  % CI include the true value?
ciWid = zeros(nTest,3);  % CI width 

% Repeat

DispProg(0, nTest, "Repeating ...");
for it=1:nTest
    [p, lm, blm, pTest, effectSize] = TestBootCompare(bTrue, eSig);
    ci = lm.coefCI;
    pErr(it,:) = log10([lm.Coefficients.pValue(2), blm.p(2), pTest]) - log10(p);
    bErr(it,:) = [lm.Coefficients.Estimate(2), blm.b(2), effectSize.mean] - bTrue;
    nGoodCI(it,:) = [prod([ci(2,1)<=bTrue, ci(2,2)>=bTrue]), ...
        prod([blm.bCI(2,1)<=bTrue, blm.bCI(2,2)>=bTrue]), ...
        prod([effectSize.ci(1)<=bTrue, effectSize.ci(2)>=bTrue]), ...
        ];
    ciWid(it,:) = [diff(ci(2,:)), diff(blm.bCI(2,:)), diff(effectSize.ci)];
    if mod(it,10) == 0
        DispProg(it, nTest);
    end
end
pErr(pErr==-Inf) = min(pErr(pErr>-Inf), [], 'all');

% Results

mtit = ["LM", "stud-resid-boot LM", "boot-test"];
xl = strings(4,3);  % [imetric im]
fig = NewFig2(2,6);
for im=1:3
    subplot(141);  histogram(pErr(:,im));  hold on;  xl(1,im) = sprintf("%.2g +/- %.2g", mean(pErr(:,im)), std(pErr(:,im)));
    subplot(142);  histogram(bErr(:,im));  hold on;  xl(2,im) = sprintf("%.2g +/- %.2g", mean(bErr(:,im)), std(bErr(:,im)));
    subplot(143);  histogram(nGoodCI(:,im));  hold on;  xl(3,im) = sprintf("%.2g +/- %.2g", mean(nGoodCI(:,im)), std(nGoodCI(:,im)));
    subplot(144);  histogram(ciWid(:,im));  hold on;  xl(4,im) = sprintf("%.2g +/- %.2g", mean(ciWid(:,im)), std(ciWid(:,im)));
end
subplot(141);  ax = gca;  ax.XLabel.String = xl(1,:);  ax.Title.String = "p error (log)";
subplot(142);  ax = gca;  ax.XLabel.String = xl(2,:);  ax.Title.String = "b error RMS";
subplot(143);  ax = gca;  ax.XLabel.String = xl(3,:);  ax.Title.String = "Num good CI";  legend(mtit, location="northwest");
subplot(144);  ax = gca;  ax.XLabel.String = xl(4,:);  ax.Title.String = "CI width";

fprintf("\n ## stud-resid-boot-LM vs boot-test: \n");
fprintf("p error (log): p=%.2g \n", BootCompareTwo(pErr(:,2), pErr(:,3), nBoot));
fprintf("b error RMS: p=%.2g \n", BootCompareTwo(bErr(:,2), bErr(:,3), nBoot));
fprintf("Num good CI: p=%.2g \n", BootCompareTwo(nGoodCI(:,2), nGoodCI(:,3), nBoot));
fprintf("CI width: p=%.2g \n", BootCompareTwo(ciWid(:,2), ciWid(:,3), nBoot));

% Findings
%   Both boot-LM and boot-test produced smaller p-values than t-test (about 0.5 order of
%   magnitude). Boot-test produced significantly better in this p error (p=0.045). No
%   difference in CI width and good CI numbers. Thus, it's better to use the boot-test.


%% Simulated data of two groups with non-normal distribution

% Simulate data

clear;
bTrue = 1.2;  % mean *ratio* of effect size 
rng('default')  % for reproducibility
yCtr = exprnd(1, 30,1);  % control group
yExp = exprnd(bTrue, 50,1);  % experimental group

% Student t-test

fig = NewFig(1,2);
subplot(121);
    histogram(yCtr);  hold on;  histogram(yExp);
    legend(["control", "experiment"])
    [~, p] = ttest2(yCtr, yExp);
    ax = gca;
    ax.XLabel.String = [sprintf("p=%.3g", p), ...
        sprintf("E(y_{ctr})=%.2f, E(y_{exp})=%.2f", mean(yCtr), mean(yExp))];
    ax.Title.String = ["Student t-test", "on original data"];
subplot(122);
    histogram(log(yCtr));  hold on;  histogram(log(yExp));
    [~, pLog] = ttest2(log(yCtr), log(yExp));
    ax = gca;
    ax.XLabel.String = sprintf("p=%.3g", pLog);
    ax.Title.String = ["Student t-test", "after log transform"];

% GLM without log transform, as a ground truth?

distr = "poisson";
linkFun = "log";
y = [yCtr; yExp];
X = [zeros(size(yCtr)); ones(size(yExp))];
glm = fitglm(X, y, Distribution=distr, Link=linkFun);
ciGLM = exp(glm.coefCI);

% Simple LM, Bootstrap LM with dummy var

nBoot = 5000;  % for min p=0.001

X = [ones(size(X)), X];  % we need the intercept column for BootLM()
[blm, lm, ciLM] = BootLM(nBoot, y, X);
bLM = sum(lm.Coefficients.Estimate)/lm.Coefficients.Estimate(1);
ciLM = sum(ciLM,1) ./ ciLM(1,:);
ciBLM = sum(blm.bCI,1) ./ blm.bCI(1,:);

% Simple LM, Bootstrap LM with dummy var *after log*

[blmLog, lmLog, ciLMlog] = BootLM(nBoot, log(y), X);
ciLMlog = exp(ciLMlog);
blmLog.b = exp(blmLog.b);  blmLog.bCI = exp(blmLog.bCI);

% Bootstrap hypothesis testing

nSided = 2;  % two-sided
effectSizeType = "meanRatio";  % meanDiff, meanRatio, medianDiff, medianRatio, or linReg
effectSizeAlpha = 0.05;
bFig = true;  % plot a related figure
[pTest, t, tBoot, effectSize] = BootCompareTwo(yCtr, yExp, nBoot, nSided, effectSizeType, effectSizeAlpha, bFig);

% Results

fprintf("\n ## meanRatio (effect size) from [true, GLM, simple-LM, bootstrap-LM, log-LM, log-boot-LM, bootstrap-test]: \n");
[bTrue, exp(glm.Coefficients.Estimate(2)), bLM, sum(blm.b)/blm.b(1), exp(lmLog.Coefficients.Estimate(2)), blmLog.b(2), effectSize.mean]
fprintf("\n ## p from [t-test, log-t-test, GLM, simple-LM, bootstrap-LM, log-LM, log-boot-LM, bootstrap-test]: \n");
[p, pLog, glm.Coefficients.pValue(2), lm.Coefficients.pValue(2), blm.p(2), lmLog.Coefficients.pValue(2), blmLog.p(2), pTest]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from GLM: \n");
[ciGLM(2,:), prod([ciGLM(2,1)<=bTrue, ciGLM(2,2)>=bTrue],2), diff(ciGLM(2,:))]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from simple-LM: \n");
[ciLM, prod([ciLM(1)<=bTrue, ciLM(2)>=bTrue],2), diff(ciLM)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from bootstrap-LM: \n");
[ciBLM, prod([ciBLM(1)<=bTrue, ciBLM(2)>=bTrue],2), diff(ciBLM)]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from log-LM: \n");
[ciLMlog(2,:), prod([ciLMlog(2,1)<=bTrue, ciLMlog(2,2)>=bTrue],2), diff(ciLMlog(2,:))]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from log-bootstrap-LM: \n");
[blmLog.bCI(2,:), prod([blmLog.bCI(2,1)<=bTrue, blmLog.bCI(2,2)>=bTrue],2), diff(blmLog.bCI(2,:))]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from GLM with %s-%s: \n", distr, linkFun);
[ciGLM(2,:), prod([ciGLM(2,1)<=bTrue, ciGLM(2,2)>=bTrue],2), diff(ciGLM(2,:))]
fprintf("\n ## [CIlow, CIhigh, CI-true, CI-width] from bootstrap-test: \n");
ci = effectSize.ci';
[ci, prod([ci(1)<=bTrue, ci(2)>=bTrue],2), diff(ci)]


%% compare GLM vs bootstrap-LM vs bootstrap-test

bTrue = 0.8;
nTest = 100;
pErr = zeros(nTest,3);  % error against t-test p in log scale
bErr = zeros(nTest,3);  
nGoodCI = zeros(nTest,3);  % CI include the true value?
ciWid = zeros(nTest,3);  % CI width 

% Repeat

DispProg(0, nTest, "Repeating ...");
for it=1:nTest
    [pLog, glm, ciGLM, blm, lm, bLM, ciLM, ciBLM, pTest, effectSize] = TestBootCompare2(bTrue);    
    pErr(it,:) = log10([glm.Coefficients.pValue(2), blm.p(2), pTest]) - log10(pLog);
    bErr(it,:) = [exp(glm.Coefficients.Estimate(2)), sum(blm.b)/blm.b(1), effectSize.mean] - bTrue;
    nGoodCI(it,:) = [prod([ciGLM(2,1)<=bTrue, ciGLM(2,2)>=bTrue]), ...
        prod([ciBLM(1)<=bTrue, ciBLM(2)>=bTrue]), ...
        prod([effectSize.ci(1)<=bTrue, effectSize.ci(2)>=bTrue]), ...
        ];
    ciWid(it,:) = [diff(ciGLM(2,:)), diff(ciBLM), diff(effectSize.ci)];
    if diff(ciBLM) < 0
        error
    end
    if mod(it,10) == 0
        DispProg(it, nTest);
    end
end
pErr(pErr==-Inf) = min(pErr(pErr>-Inf), [], 'all');

% Results

mtit = ["GLM", "stud-resid-boot LM", "boot-test"];
xl = strings(4,3);  % [imetric im]
fig = NewFig2(2,6);
for im=1:3
    subplot(141);  histogram(pErr(:,im));  hold on;  xl(1,im) = sprintf("%.2g +/- %.2g", mean(pErr(:,im)), std(pErr(:,im)));
    subplot(142);  histogram(bErr(:,im));  hold on;  xl(2,im) = sprintf("%.2g +/- %.2g", mean(bErr(:,im)), std(bErr(:,im)));
    subplot(143);  histogram(nGoodCI(:,im));  hold on;  xl(3,im) = sprintf("%.2g +/- %.2g", mean(nGoodCI(:,im)), std(nGoodCI(:,im)));
    subplot(144);  histogram(ciWid(:,im));  hold on;  xl(4,im) = sprintf("%.2g +/- %.2g", mean(ciWid(:,im)), std(ciWid(:,im)));
end
subplot(141);  ax = gca;  ax.XLabel.String = xl(1,:);  ax.Title.String = "p error (log)";
subplot(142);  ax = gca;  ax.XLabel.String = xl(2,:);  ax.Title.String = "b error RMS";
subplot(143);  ax = gca;  ax.XLabel.String = xl(3,:);  ax.Title.String = "Num good CI";  legend(mtit, location="northwest");
subplot(144);  ax = gca;  ax.XLabel.String = xl(4,:);  ax.Title.String = "CI width";

fprintf("\n ## GLM vs boot-test: \n");
fprintf("p error (log): p=%.2g \n", BootCompareTwo(pErr(:,1), pErr(:,3), nBoot));
fprintf("b error RMS: p=%.2g \n", BootCompareTwo(bErr(:,1), bErr(:,3), nBoot));
fprintf("Num good CI: p=%.2g \n", BootCompareTwo(nGoodCI(:,1), nGoodCI(:,3), nBoot));
fprintf("CI width: p=%.2g \n", BootCompareTwo(ciWid(:,1), ciWid(:,3), nBoot));

fprintf("\n ## stud-resid-boot-LM vs boot-test: \n");
fprintf("p error (log): p=%.2g \n", BootCompareTwo(pErr(:,2), pErr(:,3), nBoot));
fprintf("b error RMS: p=%.2g \n", BootCompareTwo(bErr(:,2), bErr(:,3), nBoot));
fprintf("Num good CI: p=%.2g \n", BootCompareTwo(nGoodCI(:,2), nGoodCI(:,3), nBoot));
fprintf("CI width: p=%.2g \n", BootCompareTwo(ciWid(:,2), ciWid(:,3), nBoot));

% Findings
%   When bTrue > 1, p value error and b error were not significantly different across the 
%   three methods. The resid-boot-LM produced significantly narrower CI than the other two, but 
%   produced significantly lower number of good CIs (higher chance for the CI to miss the
%   true value). 
%   When bTrue < 1, p value error was significantly lower in boot-test (good), with
%   similar good CI numbers as GLM. Although the boot-LM produced a higher good CI number
%   than boot-test, it also produced much broader CIs. 
%   Thus, the boot-test was overall better.


%% Validate with estimationstats.com

% Simulate data

clear;
bTrue = 0.2;  % mean *difference*
rng('default')  % for reproducibility
rng(23)
yCtr = exprnd(1, 30,1);  % control group
yExp = exprnd(1+bTrue, 50,1);  % experimental group

% Simple GLM with identity link

distr = "poisson";
linkFun = "identity";
y = [yCtr; yExp];
X = [zeros(size(yCtr)); ones(size(yExp))];
glm = fitglm(X, y, Distribution=distr, Link=linkFun)
ciGLM = coefCI(glm)

% Bootstrap hypothesis testing

nBoot = 5000;  % they used 5000
nSided = 2;  % two-sided
effectSizeType = "meanDiff";  % meanDiff, meanRatio, medianDiff, medianRatio, or linReg
effectSizeAlpha = 0.05;
bFig = true;  % plot a related figure
[pTest, t, tBoot, effectSize] = BootCompareTwo(yCtr, yExp, nBoot, nSided, effectSizeType, effectSizeAlpha, bFig);
effectSize
effectSize.ci'
pTest

% Findings
%   My function produced a similar or slightly wider CI than the site's result. 
%   The p value was very close too.

%% FUNCTIONS

function [p, lm, blm, pTest, effectSize] = TestBootCompare(bTrue, eSig)

    rng('shuffle')  % for varying results
    yCtr = normrnd(1, eSig, 200,1);  % control group
    yExp = normrnd(1+bTrue, eSig, 100,1);  % experimental group

    % Student t-test
    [~, p] = ttest2(yCtr, yExp);

    % Simple LM, Bootstrap LM with dummy var
    nBoot = 5000;
    y = [yCtr; yExp];
    X = [zeros(size(yCtr)); ones(size(yExp))];
    X = [ones(size(X)), X];
    bResampResid = true;  % since X is fixed.
    [blm, lm] = BootLM(nBoot, y, X, bResampResid);  % see egBootLM.m for detailed/default options

    % Bootstrap hypothesis testing

    nSided = 2;  % two-sided
    effectSizeType = "meanDiff";  % meanDiff, meanRatio, medianDiff, medianRatio, or linReg
    [pTest, t, tBoot, effectSize] = BootCompareTwo(yCtr, yExp, nBoot, nSided, effectSizeType);

end


function [pLog, glm, ciGLM, blm, lm, bLM, ciLM, ciBLM, pTest, effectSize] = ...
    TestBootCompare2(bTrue)

    rng('shuffle')  % for varying results
    yCtr = exprnd(1, 30,1);  % control group
    yExp = exprnd(bTrue, 50,1);  % experimental group

    % Student t-test *after log*
    [~, pLog] = ttest2(log(yCtr), log(yExp));

    % GLM without log transform, as a ground truth?
    distr = "poisson";
    linkFun = "log";
    y = [yCtr; yExp];
    X = [zeros(size(yCtr)); ones(size(yExp))];
    glm = fitglm(X, y, Distribution=distr, Link=linkFun);
    ciGLM = exp(glm.coefCI);

    % Simple LM, Bootstrap LM with dummy var
    nBoot = 5000;  % for min p=0.001
    X = [ones(size(X)), X];  % we need the intercept column for BootLM()
    [blm, lm, ciLM] = BootLM(nBoot, y, X);
    bLM = sum(lm.Coefficients.Estimate)/lm.Coefficients.Estimate(1);
    ciLM = sum(ciLM,1) ./ ciLM(1,:);
    ciBLM = sum(blm.bCI,1) ./ blm.bCI(1,:);

    % Simple LM, Bootstrap LM with dummy var *after log*
%     [blmLog, lmLog, ciLMlog] = BootLM(nBoot, log(y), X);
%     ciLMlog = exp(ciLMlog);
%     blmLog.b = exp(blmLog.b);  blmLog.bCI = exp(blmLog.bCI);

    % Bootstrap hypothesis testing
    nSided = 2;  % two-sided
    effectSizeType = "meanRatio";  % meanDiff, meanRatio, medianDiff, medianRatio, or linReg
    effectSizeAlpha = 0.05;
    [pTest, t, tBoot, effectSize] = BootCompareTwo(yCtr, yExp, nBoot, nSided, effectSizeType, effectSizeAlpha);

end
