%% MO-ASMO-II :: Fan2017CMOP3 -- Test Problem 3
% * TYPE-II constraint
% * Continuous, non-convex Pareto set
%
% Zhun Fan, Wenji Li, Xinye Cai, Huibiao Lin, Kaiwen Hu, Haibin Yin,
% "Difficulty Controllable and Scalable Constrained Multi-objective Test Problems,"
% 2015 International Conference on Industrial Informatics-Computing Technology, Intelligent
% Technology, Industrial Information Integration
% December 3-4, 2015, DOI: 10.1109/ICIICII.2015.105
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function Fan2017CMOP3
    % Problem setup
    problem.parameter.a = 0.51;
    problem.parameter.b = 0.5;
    problem.functions.hifi_combined_exp = @hff_combined;
    problem.functions.hifi_expensive = false;
    problem.functions.hifi_parallel = true; % if parallel pool > 1, run in parallel
    problem.functions.hifi_vectorized = true; % if no parallel pool, run in vectorized way
    problem.bound.num_x = 30;
    problem.bound.num_f = 2;
    problem.bound.xlb = zeros(1, problem.bound.num_x);
    problem.bound.xub = ones(1, problem.bound.num_x);
    problem.sampling.initial.number = 10;
    problem.surrogate.method = 'GPR';
    % Run MO-ASMO
    problem.control.casefile = mfilename('fullpath');
    result1 = runMOASMO(problem);
    % Run NSGA-II
    problem = result1.problem;
    result2 = runDO(problem, 'NSGA-II');
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [f, c, ceq] = hff_combined(x, p)
    if ((size(x, 1) ~= 1) && (size(x, 2) == 1)), x = x'; end % transpose if x is column vector
    number = size(x, 1);
    g1 = zeros(number, 1);
    g2 = zeros(number, 1);

    for idxg1 = 3:2:29
        g1 = g1 + (x(:, idxg1) - sin(0.5*pi*x(:, 1))).^2;
    end
    for idxg2 = 2:2:30
        g2 = g2 + (x(:, idxg2) - cos(0.5*pi*x(:, 1))).^2;
    end

    f1 = x(:, 1) + g1;
    f2 = 1 - x(:, 1).^2 + g2;
    f = [f1, f2];

    c1 = (g1 - p.a).*(g1 - p.b);
    c2 = (g2 - p.a).*(g2 - p.b);
    c = [c1, c2];

    ceq = [];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
