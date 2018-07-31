%% MO-ASMO-II :: runDO function
% 1. Read problem structure
% 2. Run direct optimization algorithm
% 3. Return result structure
% Usage:
%  result = runMOASMO(problem, solver)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function result = runDO(problem, solver, varargin)
    declareGlobalVariables;

    switch lower(solver)
    case 'nsga-ii'
        % Set initial populations
        if (nargin >= 3)
            initpop = varargin{1};
        else
            initpop = [];
        end
        % Combined function vs separate functions
        if (size(problem.functions.hifi_combined_exp, 1) ~= 0)
            % Combined function
            [xopt, fopt, out] = runCombinedObjNSGA2(problem, initpop);
        elseif (size(problem.functions.hifi_obj_exp, 1) ~= 0)
            % Separate function
            [xopt, fopt, out] = runSeparateObjNSGA2(problem, initpop);
        else
            error('No objective function provided');
        end
    otherwise
        error([solver, ' not supported']);
    end

    result.problem = problem;
    result.DO.xopt = xopt;
    result.DO.fopt = fopt;
    result.DO.out = out;

end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xopt, fopt, out] = runCombinedObjNSGA2(problem, initpop)
    xLast = [];
    myf = [];
    myc = [];
    myceq = [];
    fun = @objfnCombinedObj;
    cfun = @constrfnCombinedObj;
    fn.combined = problem.functions.hifi_combined_exp;
    fn.nonlcon_cheap = problem.functions.hifi_nonlcon_cheap;
    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    lb = problem.bound.xlb;
    ub = problem.bound.xub;
    num_x = problem.bound.num_x;
    
    npool = parallelPoolSize();
    opt = gaoptimset(@gamultiobj);
    if ((npool > 1) && problem.functions.hifi_parallel)
        opt.Vectorized = 'off';
        opt.UseParallel = true;
    else
        opt.UseParallel = false;
        if problem.functions.hifi_vectorized,
            opt.Vectorized = 'on';
        else
            opt.Vectorized = 'off';
        end
    end
    popsize = 10*problem.optimization.nsga2.popsize;
    if size(initpop, 1) == 0
        initpopx = [];
        initpopf = [];
    else
        [initpopx, initpopf, ~] = ndSort(initpop.x, initpop.f);
        if (size(initpopx, 1) > popsize)
            initpopx = initpopx(1:popsize, :);
            initpopf = initpopf(1:popsize, :);
        end
    end
    opt.PopulationSize = popsize;
    opt.InitialPopulation = initpopx;
    opt.InitialScores = initpopf;
    opt.ParetoFraction = problem.optimization.nsga2.paretofrac;
    opt.StallGenLimit = problem.optimization.nsga2.stallgenlimit;
    opt.PlotFcns = @gaplotpareto;

    % Call gamultiobj (NSGA-II) solver
    [xopt, fopt, exitflag, output] = gamultiobj( ...
        @(x) fun(x, problem.parameter, fn), num_x, A, b, Aeq, beq, lb, ub, ...
        @(x) cfun(x, problem.parameter, fn), opt);
    out.exitflag = exitflag;
    out.output = output;

    % ---1---------2---------3---------4---------5---------6---------7---------8---------9---------0

    function f = objfnCombinedObj(x, p, fn)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallCombinedObj(x, p, fn);
            xLast = x;
        end
        f = myf;
    end

    % ---1---------2---------3---------4---------5---------6---------7---------8---------9---------0

    function [c, ceq] = constrfnCombinedObj(x, p, fn)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallCombinedObj(x, p, fn);
            xLast = x;
        end
        c = myc;
        ceq = myceq;
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [f, c, ceq] = computeallCombinedObj(x, p, fn)
    [f, c1, ceq1] = feval(fn.combined, x, p);
    if size(fn.nonlcon_cheap, 1) == 0
        c2 = []; ceq2 = [];
    else
        [c2, ceq2] = feval(fn.nonlcon_cheap, x, p);
    end
    c = [c1, c2];
    ceq = [ceq1, ceq2];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xopt, fopt, out] = runSeparateObjNSGA2(problem, initpop)
    objfn = problem.functions.hifi_obj_exp;
    fn.nonlcon1 = problem.functions.hifi_nonlcon_exp;
    fn.nonlcon2 = problem.functions.hifi_nonlcon_cheap;
    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    lb = problem.bound.xlb;
    ub = problem.bound.xub;
    num_x = problem.bound.num_x;
    
    npool = parallelPoolSize();
    opt = gaoptimset(@gamultiobj);
    if ((npool > 1) && problem.functions.hifi_parallel)
        opt.Vectorized = 'off';
        opt.UseParallel = true;
    else
        opt.UseParallel = false;
        if problem.functions.hifi_vectorized,
            opt.Vectorized = 'on';
        else
            opt.Vectorized = 'off';
        end
    end
    popsize = problem.optimization.nsga2.popsize;
    if size(initpop, 1) == 0
        initpopx = [];
        initpopf = [];
    else
        [initpopx, initpopf, ~] = ndSort(initpop.x, initpop.f);
        if (size(initpopx, 1) > popsize)
            initpopx = initpopx(1:popsize, :);
            initpopf = initpopf(1:popsize, :);
        end
    end
    opt.PopulationSize = popsize;
    opt.InitialPopulation = initpopx;
    opt.InitialScores = initpopf;
    opt.ParetoFraction = problem.optimization.nsga2.paretofrac;
    opt.StallGenLimit = 10*problem.optimization.nsga2.stallgenlimit;
    opt.PlotFcns = @gaplotpareto;

    % Call gamultiobj (NSGA-II) solver
    if (size(fn.nonlcon1, 1) ~= 0) && (size(fn.nonlcon2, 1) ~= 0)
        [xopt, fopt, exitflag, output] = gamultiobj( ...
            @(x) objfn(x, problem.parameter), num_x, A, b, Aeq, beq, lb, ub, ...
            @(x) constrfnSeparateObj(x, problem.parameter, fn), opt);
    elseif (size(fn.nonlcon1, 1) ~= 0) && (size(fn.nonlcon2, 1) == 0)
        nonlconfun = fn.nonlcon1;
        [xopt, fopt, exitflag, output] = gamultiobj( ...
            @(x) objfn(x, problem.parameter), num_x, A, b, Aeq, beq, lb, ub, ...
            @(x) nonlconfun(x, problem.parameter), opt);
    elseif (size(fn.nonlcon1, 1) == 0) && (size(fn.nonlcon2, 1) ~= 0)
        nonlconfun = fn.nonlcon2;
        [xopt, fopt, exitflag, output] = gamultiobj( ...
            @(x) objfn(x, problem.parameter), num_x, A, b, Aeq, beq, lb, ub, ...
            @(x) nonlconfun(x, problem.parameter), opt);
    elseif (size(fn.nonlcon1, 1) == 0) && (size(fn.nonlcon2, 1) == 0)
        [xopt, fopt, exitflag, output] = gamultiobj( ...
            @(x) objfn(x, problem.parameter), num_x, A, b, Aeq, beq, lb, ub, [], opt);
    end
    out.exitflag = exitflag;
    out.output = output;
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [c, ceq] = constrfnSeparateObj(x, p, fn)
    [c1, ceq1] = feval(fn.nonlcon1, x, p);
    [c2, ceq2] = feval(fn.nonlcon2, x, p);
    c = [c1, c2];
    ceq = [ceq1, ceq2];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
