%% MO-ASMO-II :: trainInvalidRegion function
% 1. Train model for regions with invalid inputs
% Usage:
%  [irmodel, t_elapsed] = trainInvalidRegion(problem, xinvalid)
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [irmodel, t_elapsed] = trainInvalidRegion(problem, xinvalid)
    declareGlobalVariables;
    
    t_elapsed = 0;
    tic;
    
    % Use invalid region model?
    if (problem.invalidregion.use == false)
        irmodel = [];
        return;
    end

    if verbose, disp('Train invalid region models...'); end

    % Scale input X (xinvalid)
    scale_x = problem.surrogate.scale;
    xlb = problem.bound.xlb;
    xub = problem.bound.xub;
    if scale_x
        xinvalid = varScale(xinvalid, xlb, xub, 'scale');
    end

    p.C = problem.invalidregion.C;
    p.q = problem.invalidregion.q;
    p.epsilon = problem.invalidregion.epsilon;
    p.x = xinvalid;
    number = size(p.x, 1);
    if (number < 5)
        irmodel = [];
        return;
    end

    d0 = ones(number, 1);
    lb = zeros(size(d0));
    ub = p.C*ones(size(d0));
    Aeq = ones(1, number);
    beq = 1;
    opt = optimoptions('fmincon', 'Algorithm', 'sqp', 'ConstraintTolerance', 1e-6, ...
        'Display', 'none', 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, ...
        'OptimalityTolerance', 1e-3, 'ScaleProblem', true, 'UseParallel', false, ...
        'StepTolerance', 1e-4);
    npool = parallelPoolSize();
    if (npool > 1), opt.UseParallel = true; end
    [dopt, ~] = fmincon(@(d) obj(d, p), d0, [], [], Aeq, beq, lb, ub, [], opt);

    irmodel.dopt = dopt;
    irmodel.parameter = p;
    irmodel.scale.scale_x = scale_x;
    irmodel.scale.xlb = xlb;
    irmodel.scale.xub = xub;
    
    t_elapsed = toc;
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function f = obj(d, p)
    x = p.x;
    npx = size(x, 1);
    f1 = 0;
    for m = 1:npx
        f1 = f1 + d(m)*exp(0);
    end
    f2 = 0;
    for n = 1:npx
        for m = 1:npx
            f2 = f2 + d(m)*d(n)*exp(-p.q*sum((x(m, :) - x(n, :)).^2));
        end
    end
    f = -(f1 - f2);
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
