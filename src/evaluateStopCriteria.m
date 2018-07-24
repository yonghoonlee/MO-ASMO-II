%% MO-ASMO-II :: evaluateStopCriteria function
% 1. Check if any or all stopping criteria are met and return true/false.
% Usage:
%  varargout = evalStopCriteria(varargin)
% Arguments:
%  {problem, k, stopcount, EDvec, resHVpred, resHVRpred, resHVhff, resHVRhff}
% Returns:
%  {stoploop, stopcount}
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [stoploop, stopcount] = evaluateStopCriteria(problem, k, stopcount, EDvec, ...
            resHVpred, resHVRpred, resHVhff, resHVRhff)
    % Default return values are false
    stoploop = false;
    % Get parameters
    maxiter = problem.stop.maxiter;
    method = problem.stop.residual.satisfaction_method;
    continuous = problem.stop.residual.satisfaction_continuous;
    ED_max = problem.stop.residual.ED_max;
    ED_avg = problem.stop.residual.ED_avg;
    HV_size = problem.stop.residual.HV_size;
    HV_ratio = problem.stop.residual.HV_ratio;
    HV_data = problem.stop.residual.HV_data;
    % Default parameters for stop criteria are false
    criteria_ED_max = false;
    criteria_ED_avg = false;
    criteria_HV = false;
    criteria_HVR = false;

    % 1. Check if maxiter reached
    if k >= maxiter
        stoploop = true;
        return;
    end

    % 2. Check if Euclidean distance criterion is met
    if (max(EDvec) <= ED_max), criteria_ED_max = true; end
    if (mean(EDvec) <= ED_avg), criteria_ED_avg = true; end

    % 3. Check if HV residual criteria are met
    switch lower(HV_data)
    case 'highfidelity'
        if (resHVhff <= HV_size), criteria_HV = true; end
        if (resHVRhff <= HV_ratio), criteria_HVR = true; end
    case 'predicted'
        if (resHVpred <= HV_size), criteria_HV = true; end
        if (resHVRpred <= HV_ratio), criteria_HVR = true; end
    case 'both'
        if ((resHVhff <= HV_size) && (resHVpred <= HV_size)), criteria_HV = true; end
        if ((resHVRhff <= HV_ratio) && (resHVRpred <= HV_ratio)), criteria_HVR = true; end
    otherwise
        error([HV_data, ' option not supported']);
    end

    % 4. Evaluate stopping condition
    switch lower(method)
    case 'and'
        if (criteria_ED_max && criteria_ED_avg && criteria_HV && criteria_HVR)
            stopcount = stopcount + 1;
        else
            stopcount = 0;
        end
    case 'or'
        if (criteria_ED_max || criteria_ED_avg || criteria_HV || criteria_HVR)
            stopcount = stopcount + 1;
        else
            stopcount = 0;
        end
    otherwise
        error([method, ' option not supported']);
    end
    if (stopcount >= continuous), stoploop = true; end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
