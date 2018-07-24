%% MO-ASMO-II :: approxNDHV function
% 1. Evaluate hypervolume and hypervolume ratio convergence metrices with given Pareto set
% Usage:
%  [HV, HVR] = approxNDHV(problem, fPareto)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [HV, HVR] = approxNDHV(problem, fPareto)
    n = problem.stop.residual.HV_eval;
    [number, num_f] = size(fPareto);

    uf = min(fPareto, [], 1); % utopia
    auf = max(fPareto, [], 1); % antiutopia

    samples = uf + (auf - uf).*rand(n, num_f);
    dominated = 0;

    for idx1 = 1:number
        idx2 = (sum((fPareto(idx1, :) <= samples), 2) == num_f);
        dominated = dominated + sum(idx2);
        samples(idx2, :) = [];
    end

    hvtotal = prod(auf - uf);
    HVR = 1 - dominated/n;
    HV = HVR*hvtotal;
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
