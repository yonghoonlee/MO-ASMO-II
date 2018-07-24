%% MO-ASMO-II :: errEuclideanDistance function
% 1. Evaluate Euclidean distance error between high fidelity and predicted results of val. points
% Usage:
%  varargout = errEuclideanDistance(varargin)
% Arguments:
%  {problem, hffF, surF, prevEDarr, fmodel}
% Returns:
%  {EDvec, EDavg}
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [EDvec, EDavg] = errEuclideanDistance(problem, hffF, surF, prevEDarr, fmodel)
    % If no surrogate exists
    if (size(fmodel, 1) == 0)
        EDvec = [];
        EDavg = [];
        return;
    end    

    if ((size(hffF, 1) == 0) && (size(prevEDarr, 1) == 0))  % k = 1, no data
        EDvec = 1;
        EDavg = 1;
    elseif (size(hffF, 1) == 0) % k > 1, no data
        maxError = 0;
        for idx = 1:size(prevEDarr, 1)
            data = prevEDarr{idx, 1};
            maxError = max(maxError, max(data));
        end
        EDvec = maxError;
        EDavg = maxError;
    else % data exists
        % Scale variable
        scale_f = fmodel.scale.scale_f;
        if scale_f
            flb = fmodel.scale.flb;
            fub = fmodel.scale.fub;
            hffF = varScale(hffF, flb, fub, 'scale');
            surF = varScale(surF, flb, fub, 'scale');
        end
        % Compute ED
        EDvec = sqrt(sum((hffF - surF).^2, 2));
        EDavg = mean(EDvec);
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
