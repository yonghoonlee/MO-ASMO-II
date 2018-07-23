%% MO-ASMO-II :: samplingValidation function
% 1. Sample validation points from predicted Pareto set
% Usage:
%  x = samplingValidation(problem, k)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [valX, valF, valC, valCEQ] = samplingValidation(problem, parX, parF, parC, parCEQ)
    declareGlobalVariables;
    if verbose, disp('Sampling for validation...'); end

    num_x = problem.bound.num_x;
    num_f = problem.bound.num_f;
    xnumber = size(parX, 1);
    Vnumber = problem.sampling.validation.number;
    Vmethod = problem.sampling.validation.method;
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
