%% MO-ASMO-II :: samplingUpdate function
% 1. Generate update samples for updating surrogate model
% Usage:
%  xt = samplingUpdate(problem, k)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function xt = samplingUpdate(problem, k, poolXvalid, irmodel, parX)
    declareGlobalVariables;

    exploit_method = problem.sampling.update.exploit.method;
    exploit_number = problem.sampling.update.exploit.number;
    explore_method = problem.sampling.update.explore.method;
    explore_number = problem.sampling.update.explore.number;
    dimension = problem.bound.num_x;
    xlb = problem.bound.xlb;
    xub = problem.bound.xub;

    % Exploration
    ex = -1;
    while(ex < 0)
        ex = 2;
        if verbose, disp('Exploration sampling...'); end

        switch lower(explore_method)
        case 'lhs'
            % random permutation to create Latin Hypercube
            xt1 = zeros(explore_number, dimension);
            for idx = 1:dimension
                xt1(:,idx) = randperm(explore_number)';
            end
            % randomize within Latin Hypercube
            xt1 = (xt1 - 0.5)/explore_number;
            rn = rand(size(xt1));
            rn = (rn - 0.5)/explore_number;
            xt1 = xt1 + rn;
        case 'random'
            xt1 = rand(explore_number, dimension);
        otherwise
            error([explore_method,' not supported']);
        end

        % Unscale
        xt1 = varScale(xt1, xlb, xub, 'unscale');

        % Adjust samples to (1) avoid existing points
        %                   (2) comply linear constraints and cheap nonlinear constraints
        npool = parallelPoolSize();
        flghist = zeros(size(xt1, 1), 1);
        if ((problem.functions.hifi_parallel == true) && (npool > 1)) % Parallel
            currentpool = gcp('nocreate');
            for idx = 1:size(xt1, 1)
                fevFuture(idx) = parfeval(currentpool, fmincon_explore, 2, ...
                    problem, xt1(idx, :), poolXvalid, irmodel);
            end
            for idx = 1:size(xt1, 1)
                [completedIdx, value, flg] = fetchNext(fevFuture);
                xt1(completedIdx, :) = reshape(value, 1, numel(value));
                ex = min(ex, flg);
                flghist(completedIdx, 1) = flg;
            end
        else % Serial
            for idx = 1:size(xt1, 1)
                [value, flg] = fmincon_explore(problem, xt(idx, :), poolXvalid, irmodel);
                xt1(idx, :) = reshape(value, 1, numel(value));
                ex = min(ex, flg);
                flghist(idx, 1) = flg;
            end
        end

        if (sum(flghist<0) <= 0.2*size(xt1, 1))
            ex = 0;
            xt1 = xt1((flghist >= 0), :);
        end
    end

    % Check number of samples generated from exploration sampling.
    % If undersampled, sample more during exploitation sampling.
    if (size(xt1, 1) < explore_number)
        exploit_number = exploit_number + (explore_number - size(xt1, 1));
    end

    % Exploitation
    

    % Combine
    xt = [xt1; xt2];
    xt = unique(xt, 'rows');

end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xout, flg] = fmincon_explore(problem, xin, xprev, irmodel)
    % Adjust samples to (1) avoid existing points
    %                   (2) comply linear constraints and cheap nonlinear constraints
    x0 = xin;
    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    xlb = problem.bound.xlb;
    xub = problem.bound.xub;
    w1 = problem.sampling.update.explore.w1_distance;
    w2 = problem.sampling.update.explore.w2_disperse;
    opt = optimoptions('fmincon', 'Algorithm', 'sqp', 'ConstraintTolerance', 1e-6, ...
        'Display', 'none', 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, ...
        'OptimalityTolerance', 1e-3, 'ScaleProblem', true, 'UseParallel', false, ...
        'StepTolerance', 1e-12);
    hifi_nonlcon_cheap = problem.functions.hifi_nonlcon_cheap;
    [xout, ~, flg] = fmincon(@(x) explore_obj(x, x0, xprev, w1, w2), xin, A, b, Aeq, beq, ...
        xlb, xub, @(x) explore_nonlcon(hifi_nonlcon_cheap, x, problem.parameter, irmodel), opt);
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function f = explore_obj(x, x0, xprev, w1, w2)
    f1 = sum((x - x0).^2);
    f2 = 0;
    if size(xprev, 1) ~= 0
        f2 = sum(-log(max(sum((xm - xprev).^2, 2), 1e-12))) / size(xprev, 1);
    end
    f = w1*f1 + w2*f2;
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [c, ceq] = explore_nonlcon(hifi_nonlcon_cheap, x, p, irmodel)
    c1 = [];
    ceq = [];
    if size(hifi_nonlcon_cheap, 1) ~= 0
        [c1, ceq] = hifi_nonlcon_cheap(x, p);
    end
    c2 = invalidRegionEval(x, irmodel);
    c = [c1, c2];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
