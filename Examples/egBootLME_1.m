% Clustered data with normally distributed residuals
%   We can use fitlme in this case, but let's confirm BootLME produces similar results.

clear;  close all;  clc;
comInit;


%% Create fake data to simulate how Math score is affected by IQ and height, when measured 
%   school by school, and how the relation differs between states.

% set parameters
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
head(tb)


%% See how LM fails

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

% check if the residuals are normally distributed
if ~IsNormDist(lm3.Residuals.Raw)
    warning("Residuals are not normally distributed.");
end
% turns out not, likely because of the school-specific random effect

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


%% See how LME succeeds

% fit with all predictor variables
lme = fitlme(tb, 'score ~ iq + height + state + (1|school)')

% remove outliers. fitlme does not provide Diagnostics property
[~, bOut] = rmoutliers(lme.residuals);
ioOut = find(bOut);  % indices of observation of outliers
NewFig2(2,2);
line(lme.fitted, lme.residuals, Marker='x', LineStyle='none');
if sum(bOut) > 0
    x = lme.fitted;  y = lme.residuals;
    line(x(bOut), y(bOut), Marker='o', Color='r', LineStyle='none')
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

% check if the residuals are normally distributed
if ~IsNormDist(lme5.residuals)
    warning("Residuals are not normally distributed.");
end
% turns out yes

% compare the true vs estimated 
tbRes = table;
tbRes.predictor = ["State B", "State C", "IQ"]';
tbRes.true = [bState(2:end)'; bIq];
tbRes.estimate = lme5.Coefficients.Estimate(2:end);
ci = lme5.coefCI;
tbRes.ciLow = ci(2:end,1);
tbRes.ciHigh = ci(2:end,2);
tbRes.p = lme5.Coefficients.pValue(2:end);
tbResLME = tbRes
% now reveals effect of both states are significant, although the CI does not include the
% true value in state C.


%% See how BootLME works

% below are inputs to BootLME
% tb : data table to be inputted to fitlme()
% md : formula to be inputted to fitlme()
% varNames : 3-element string array:
%               varNames(1) : y variable name
%               varNames(2) : group variable name. When groups are determined by
%                             more than one variable, make a column of a unique
%                             group identifier and indicate it here. 
%               varNames(3) : cluster/subject variable name
%               e.g,. varNames = ["y", "grpId", "subjId"]
% oResamp : How to resample.    
%       'resid'     = [default] Resample all residuals, randomly across groups and subjects.
%                     Residulas can move to different groups.
%                     Good in most cases, including longitudinal data where x or time has 
%                     meaning (should be maintained)
%       'residGrp'  = Resample residuals within each group, randomly across clusters.
%                     Keeps the group structure. 
%                     Not recommended unless residuals are extremely different between groups.
%       'obsGrp'    = Resample observations (original y and x values) within each group. 
%                     Keeps the group structure. Observation number per cluster can vary.
%                     Good when measurements are randomly repeated
%       'obsClu'    = Resample observations (original y and x values) within each *cluster/subject*. 
%                     Keeps the group-cluster structure. Observation number per cluster unchanged.
%       'cluGrp'    = Resample clusters/subjects (keeping their original y and X) within each group.
%                     Keep the group structure.
%                     Good for longitudinal data with many subjects per group, where x or
%                     time within each subject should be maintained.
%               Whether to resample residuals only (true) or resample y & X? 
%               Resample only residuals when X is fixed (to retain the structure of X data). For example, X is
%               fixed in *controlled* experiment (e.g., you measure dose-response curve from a set of
%               evenly-distributed dose values). In contrast, X is not fixed in
%               *survey* experiment (e.g., you analyze relation between height (x) and 
%               weight (y) of people: when you re-recruit (resample) people, the height (x) will vary. 
% ciType : only "norm" or "per" for LME
% ciAlpha : 0.05 (default)
% weight : The option Weights of bootstrp

% use the best model found above with fitlme
nBoot = 2000;  % for p value resolution of 0.001
md = "score ~ iq + state + (1|school)";
varNames = ["score", "state", "school"];  % y, group, cluster

% bootstrap coeffiicents by resampling residuals
[blme, lme, fig] = BootLME(nBoot, tb, md, varNames, ...  % required inputs
    bFig=true);  % options for bootstrap
tbRes.estimate = blme.b(2:end);  tbRes.ciLow = blme.ci(2:end,1);  tbRes.ciHigh = blme.ci(2:end,2);  tbRes.p = blme.p(2:end);
tbRes
% all the bootstrapped coefficients are normally distributed (central limit theorem)
% very similar to LME

% test if resampling observations within each group improves the result
[blme1, ~, fig1] = BootLME(nBoot, tb, md, varNames, ...  % required inputs
    oResamp='obsGrp', bFig=true);  % options for bootstrap
tbRes1 = tbRes;
tbRes1.estimate = blme1.b(2:end);  tbRes1.ciLow = blme1.ci(2:end,1);  tbRes1.ciHigh = blme1.ci(2:end,2);  tbRes1.p = blme1.p(2:end);
tbRes1
% similar to the default BootLME above, but with broader CI (state B CI now includes zero)
    
% test if resampling observations within each cluster improves the result
[blme2, ~, fig2] = BootLME(nBoot, tb, md, varNames, ...  % required inputs
    oResamp='obsClu', bFig=true);  % options for bootstrap
tbRes2 = tbRes;
tbRes2.estimate = blme2.b(2:end);  tbRes2.ciLow = blme2.ci(2:end,1);  tbRes2.ciHigh = blme2.ci(2:end,2);  tbRes2.p = blme2.p(2:end);
tbRes2
% this produced the narrowest CI, making the error of (true - ciHigh) bigger.
    
% test if resampling schools (keeping score values within each school unchanged) imrpves
% the result
[blme3, ~, fig3] = BootLME(nBoot, tb, md, varNames, ...  % required inputs
    oResamp='cluGrp', bFig=true);  % options for bootstrap
tbRes3 = tbRes;
tbRes3.estimate = blme3.b(2:end);  tbRes3.ciLow = blme3.ci(2:end,1);  tbRes3.ciHigh = blme3.ci(2:end,2);  tbRes3.p = blme3.p(2:end);
tbRes3
% this made both state effects insignificant (wrong).

% repeat the best again with different fitlme options
[blme4, lme4, fig4] = BootLME(nBoot, tb, md, varNames, ...  % required inputs
    bFig=true, ...  % options for bootstrap
    CovariancePattern='Diagonal', FitMethod='REML');  % options to be passed to fitlme()
tbRes4 = tbRes;
tbRes4.estimate = blme4.b(2:end);  tbRes4.ciLow = blme4.ci(2:end,1);  tbRes4.ciHigh = blme4.ci(2:end,2);  tbRes4.p = blme4.p(2:end);
tbRes4
% almost no change from tbRes


