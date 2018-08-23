%% MO-ASMO-II :: scriptHVResidualParetoHff script
% 1. Compute residuals of hypervolume of true Pareto set (high-fidelity function results)
% Usage:
%  scriptHVResidualParetoHff
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

if k == 1
    c37_minHVhff = c36_HVhff;
    c37_maxHVhff = c36_HVhff;
    c37_resHVhff = 1;
    c37_resHVhffHistory = 1;
    c37_minHVRhff = c36_HVRhff;
    c37_maxHVRhff = c36_HVRhff;
    c37_resHVRhff = 1;
    c37_resHVRhffHistory = 1;
else
    c37_minHVhff = min(c37_minHVhff, c36_HVhff);
    c37_maxHVhff = max(c37_maxHVhff, c36_HVhff);
    c37_resHVhff = abs(c36_HVhff - c36_HVhffHistory(end-1,1)) ...
        / abs(c37_maxHVhff - c37_minHVhff);
    c37_resHVhffHistory = [c37_resHVhffHistory; c37_resHVhff];
    c37_minHVRhff = min(c37_minHVRhff, c36_HVRhff);
    c37_maxHVRhff = max(c37_maxHVRhff, c36_HVRhff);
    c37_resHVRhff = abs(c36_HVRhff - c36_HVRhffHistory(end-1,1)) ...
        / abs(c37_maxHVRhff - c37_minHVRhff);
    c37_resHVRhffHistory = [c37_resHVRhffHistory; c37_resHVRhff];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0