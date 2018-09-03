%% MO-ASMO-II :: scriptHVResidualParetoPred script
% 1. Compute hypervolume of predicted Pareto set and its residual
% Usage:
%  scriptHVResidualParetoPred
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

% Hypervolume (HV) and hypervolume ratio (HVR) of predicted Pareto set
[c34_HVpred, c34_HVRpred] = approxNDHV(problem, c19_parSurF_valid);
if k == 1
    c34_HVpredHistory = c34_HVpred;
    c34_HVRpredHistory = c34_HVRpred;
else
    c34_HVpredHistory = [c34_HVpredHistory; c34_HVpred];
    c34_HVRpredHistory = [c34_HVRpredHistory; c34_HVRpred];
end

% Normalized residuals of HV and HVR of predicted Pareto set
if k == 1
    c35_minHVpred = c34_HVpred;
    c35_maxHVpred = c34_HVpred;
    c35_resHVpred = 1;
    c35_resHVpredHistory = 1;
    c35_minHVRpred = c34_HVRpred;
    c35_maxHVRpred = c34_HVRpred;
    c35_resHVRpred = 1;
    c35_resHVRpredHistory = 1;
else
    c35_minHVpred = min(c35_minHVpred, c34_HVpred);
    c35_maxHVpred = max(c35_maxHVpred, c34_HVpred);
    c35_resHVpred = abs(c34_HVpred - c34_HVpredHistory(end-1,1)) ...
        / abs(c35_maxHVpred - c35_minHVpred);
    c35_resHVpredHistory = [c35_resHVpredHistory; c35_resHVpred];
    c35_minHVRpred = min(c35_minHVRpred, c34_HVRpred);
    c35_maxHVRpred = max(c35_maxHVRpred, c34_HVRpred);
    c35_resHVRpred = abs(c34_HVRpred - c34_HVRpredHistory(end-1,1)) ...
        / abs(c35_maxHVRpred - c35_minHVRpred);
    c35_resHVRpredHistory = [c35_resHVRpredHistory; c35_resHVRpred];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
