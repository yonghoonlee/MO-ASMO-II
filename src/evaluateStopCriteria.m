%% MO-ASMO-II :: evaluateStopCriteria function
% 1. Check if any or all stopping criteria are met and return true/false.
% Usage:
%  varargout = evalStopCriteria(varargin)
% Arguments:
%  {problem, k, stopcount, EDvec, resHVpred, resHVRpred, resHVhff, resHVRhff}
% Returns:
%  {stoploop, stopcount}
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [stoploop, stopcount] = evaluateStopCriteria(problem, k, stopcount, EDvec, ...
            resHVpred, resHVRpred, resHVhff, resHVRhff)
    declareGlobalVariables;
    
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
    HV_type = problem.stop.residual.HV_type;
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
    val_ED_max = max(EDvec);
    val_ED_avg = mean(EDvec);
    if (val_ED_max <= ED_max), criteria_ED_max = true; end
    if (val_ED_avg <= ED_avg), criteria_ED_avg = true; end

    % 3. Check if HV residual criteria are met
    switch lower(HV_data)
    case 'highfidelity'
        if (resHVhff <= HV_size), criteria_HV = true; end
        if (resHVRhff <= HV_ratio), criteria_HVR = true; end
        switch lower(HV_type)
        case 'absolute'
            criteria_HVR = true; % do not check ratio value
        case 'ratio'
            criteria_HV = true; % do not check absolute value
        case 'both'
            % do nothing
        end
    case 'predicted'
        if (resHVpred <= HV_size), criteria_HV = true; end
        if (resHVRpred <= HV_ratio), criteria_HVR = true; end
        switch lower(HV_type)
        case 'absolute'
            criteria_HVR = true; % do not check ratio value
        case 'ratio'
            criteria_HV = true; % do not check absolute value
        case 'both'
            % do nothing
        end
    case 'both'
        if ((resHVhff <= HV_size) && (resHVpred <= HV_size)), criteria_HV = true; end
        if ((resHVRhff <= HV_ratio) && (resHVRpred <= HV_ratio)), criteria_HVR = true; end
        switch lower(HV_type)
        case 'absolute'
            criteria_HVR = true; % do not check ratio value
        case 'ratio'
            criteria_HV = true; % do not check absolute value
        case 'both'
            % do nothing
        end
    otherwise
        error([HV_data, ' option not supported']);
    end

    % 4. Evaluate stopping condition
    switch lower(method)
    case 'and'
        if (criteria_ED_max && criteria_ED_avg && (criteria_HV && criteria_HVR))
            stopcount = stopcount + 1;
            if verbose
                disp(['Stopping condition satisfied for ', num2str(stopcount), ' iterations']);
            end
        else
            stopcount = 0;
        end
    case 'or'
        if (criteria_ED_max || criteria_ED_avg || (criteria_HV && criteria_HVR))
            stopcount = stopcount + 1;
            if verbose
                disp(['Stopping condition satisfied for ', num2str(stopcount), ' iterations']);
            end
        else
            stopcount = 0;
        end
    otherwise
        error([method, ' option not supported']);
    end
    
    % 5. Display values
    if verbose
        disp([sprintf('%12.4e',val_ED_max), sprintf('%12.4e',val_ED_avg), ...
            sprintf('%12.4e',resHVpred), sprintf('%12.4e',resHVRpred), ...
            sprintf('%12.4e',resHVhff), sprintf('%12.4e',resHVRhff)]);
    end
    
    % 6. Stop if termination condition is met
    if (stopcount >= continuous)
        stoploop = true;
        if verbose, disp('Terminating MO-ASMO'); end
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
