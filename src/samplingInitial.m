%% MO-ASMO-II :: samplingInitial function
% 1. 
% Usage:
%  x = samplingInitial(problem, k)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function x = samplingInitial(problem)
    declareGlobalVariables;

    method = problem.sampling.initial.method;
    number = problem.sampling.initial.number;
    force_number = problem.sampling.initial.force_number;
    dimension = problem.bound.num_x;
    xlb = problem.bound.xlb;
    xub = problem.bound.xub;

    ex = -1;
    while(ex < 0)
        ex = 2;
        if verbose, disp('Sampling...'); end

        switch lower(method)
        case 'lhs'
            % random permutation to create Latin Hypercube
            xt = zeros(number, dimension);
            for idx = 1:dimension
                xt(:,idx) = randperm(number)';
            end
            % randomize within Latin Hypercube
            xt = (xt - 0.5)/number;
            rn = rand(size(xt));
            rn = (rn - 0.5)/number;
            xt = xt + rn;
        case 'random'
            xt = rand(number, dimension);
        case 'cci'
            if (dimension < 2), error(['Variables in x need to be 2 or larger for ', method]); end
            if (dimension > 10), error(['Too many variables in x for ', method]); end
            xt = (ccdesign(dimension, 'Type', 'Inscribed') + 1.0)./2;
            xt = unique(xt, 'rows');
        case 'bbd'
            if (dimension < 3), error(['Variables in x need to be 3 or larger for ', method]); end
            if (dimension > 10), error(['Too many variables in x for ', method]); end
            xt = (bbdesign(dimension) + 1.0)./2;
            xt = unique(xt, 'rows');
        otherwise
            error([method,' not supported']);
        end

        % Descale
        xt = samplingScale(xt, xlb, xub, 'descale');
        
        % Adjust samples to comply linear constraints and cheap nonlinear constraints
        npool = parallelPoolSize();
        if ((problem.functions.hifi_parallel == true) && (npool > 1)) % Parallel
            currentpool = gcp('nocreate');
            for idx = 1:size(xt,1)
                fevFuture(idx) = parfeval(currentpool, hf_fmincon, 2, xt(idx,:), problem);
            end
            for idx = 1:size(xt,1)
                [completedIdx, value, flg] = fetchNext(fevFuture);
                xt(completedIdx,:) = reshape(value, 1, numel(value));
                ex = min(ex, flg);
            end
        else % Serial
            for idx = 1:size(xt,1)
                [value, flg] = hf_fmincon(xt(idx,:), problem);
                xt(idx,:) = reshape(value, 1, numel(value));
                ex = min(ex, flg);
            end
        end

        % Force number of samples
        if (force_number == true) && (number < size(xt,1))
            xt = datasample(xt, number, 'Replace', false);
        end

    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xout, flg] = hf_fmincon(xin, problem)
    x0 = xin;
    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    xlb = problem.bound.xlb;
    xub = problem.bound.xub;
    opt = optimoptions('fmincon', 'Algorithm', 'sqp', 'ConstraintTolerance', 1e-6, ...
        'Display', 'none', 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, ...
        'OptimalityTolerance', 1e-3, 'ScaleProblem', true, 'UseParallel', false, ...
        'StepTolerance', 1e-4);
    hifi_nonlcon_cheap = problem.functions.hifi_nonlcon_cheap;
    [xout, ~, flg] = fmincon(@(x) hf_fmincon_obj(x,x0), xin, A, b, Aeq, beq, xlb, xub, ...
        @(x) hifi_nonlcon_cheap(x,p), opt);
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function f = hf_fmincon_obj(x, x0)
    f = sum((x - x0).^2);
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
