%% MO-ASMO-II :: scriptHVResidualParetoHff script
% 1. Compute hypervolume of true Pareto set (high-fidelity function results) and its residual
% Usage:
%  scriptHVResidualParetoHff
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

% Hypervolume (HV) and hypervolume ratio (HVR) of Pareto set of high fidelity results
HffC = max(c08_poolHffC_valid,[],2);
HffCEQ = max(abs(c08_poolHffCEQ_valid),[],2);
iHff = ones(size(c07_poolX_valid,1),1);
iHff(HffC>problem.control.tolC) = 0;
iHff(HffCEQ>problem.control.tolCEQ) = 0;
iHff = logical(iHff);
iHff = enforceIndexLincon(problem, iHff, c07_poolX_valid);
[~,ndFHF,ndiHF] = ndSort(c07_poolX_valid(iHff == 1,:), c08_poolHffF_valid(iHff == 1,:));
[c36_HVhff, c36_HVRhff] = approxNDHV(problem, ndFHF(ndiHF == 1,:));
if k == 1
    c36_HVhffHistory = c36_HVhff;
    c36_HVRhffHistory = c36_HVRhff;
else
    c36_HVhffHistory = [c36_HVhffHistory; c36_HVhff];
    c36_HVRhffHistory = [c36_HVRhffHistory; c36_HVRhff];
end

% Normalized residuals of HV and HVR of Pareto set of high fidelity results
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
