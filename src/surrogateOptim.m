%% MO-ASMO-II :: surrogateOptim function
% 1. Using surrogate models for f, c, ceq, run surrogate-based optimization
% Usage:
%  varargout = surrogateOptim(varargin)
% Arguments:
%  {problem, surrF, surrC, surrCEQ, irmodel, startpts}
% Returns:
%  {xP, fP, cP, ceqP, out, t_elapsed}
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xP, fP, cP, ceqP, out, t_elapsed] = surrogateOptim( ...
        problem, surrF, surrC, surrCEQ, irmodel, startpts)
    declareGlobalVariables;
    
    t_elapsed = 0;
    tic;
    
    startptX = startpts.x;
    startptF = startpts.f;
    
    solver = problem.optimization.solver;
    if verbose, disp('Solving multiobjective optimization using surrogate models'); end

    num_x = problem.bound.num_x;
    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    lb = reshape(problem.bound.xlb, 1, num_x);
    ub = reshape(problem.bound.xub, 1, num_x);

    xP = [];
    fP = [];
    cP = [];
    ceqP = [];
    out = [];

    switch lower(solver)
    case 'nsga-ii'
        popsize = problem.optimization.nsga2.popsize;
        paretofrac = problem.optimization.nsga2.paretofrac;
        stallgenlimit = problem.optimization.nsga2.stallgenlimit;
        maxgen = problem.optimization.nsga2.maxgen;
        functiontolerance = problem.optimization.nsga2.functiontolerance;

        opt = gaoptimset(@gamultiobj);
        opt.Vectorized = 'on';
        opt.UseParallel = false;
        opt.PopulationSize = popsize;
        opt.ParetoFraction = paretofrac;
        opt.StallGenLimit = stallgenlimit;
        opt.InitialPopulation = startptX;
%        opt.InitialScores = startptF;
        opt.PlotFcns = @gaplotpareto;
        opt.Generations = maxgen;
        opt.TolFun = functiontolerance;

        [xopt, fopt, exitflag, output] = gamultiobj( ...
            @(x) surrogateFeval(x, surrF), num_x, A, b, Aeq, beq, lb, ub, ...
            @(x) surrogateCeval(x, surrC, surrCEQ, irmodel), opt);
        
        out.exitflag = exitflag;
        out.output = output;
    case 'epsilon-constraints'
        nMS = problem.optimization.EC.numMultiStart;
        MS = MultiStart;
        npool = parallelPoolSize();
        if (npool > 1), MS.UseParallel = true; end
        if verbose, MS.Display = 'iter'; else MS.Display = 'off'; end
        ECprob = createOptimProblem('fmincon');
        ECprob.Aineq = problem.lincon.A;
        ECprob.bineq = problem.lincon.b;
        ECprob.Aeq = problem.lincon.Aeq;
        ECprob.beq = problem.lincon.beq;
        ECprob.lb = problem.bound.xlb;
        ECprob.ub = problem.bound.xub;
        ECprob.options.Algorithm = problem.optimization.fmincon.solver;
        ECprob.options.MaxFunctionEvaluations = Inf;
        
        % Find anchor points
        mx = problem.bound.num_x;
        mf = problem.bound.num_f;
        xopt_anchor = zeros(mf, mx);
        fopt_anchor = zeros(mf, mf);
        eopt_anchor = zeros(mf, 1);
        for idx = 1:mf
            [~, idxsrt] = sortrows(startptF, idx);
            xprevtmp = startptX(idxsrt, :);
            if size(xprevtmp, 1) == 0
                xprevtmp = (problem.bound.xlb + problem.bound.xub)./2;
            end
            ECprob.x0 = xprevtmp(1, :);
            if nMS < size(xprevtmp, 1)
                x0arr = xprevtmp(1:nMS, :);
            else
                x0arr = xprevtmp;
            end
            tpoints = CustomStartPointSet(x0arr);
            epscon = [];
            [xoptAP, foptAP, ooptAP] = runECs(problem, MS, ECprob, tpoints, ...
                idx, epscon, surrF, surrC, surrCEQ, irmodel);
            xopt_anchor(idx, :) = xoptAP;
            fopt_anchor(idx, :) = foptAP;
            eopt_anchor(idx, 1) = ooptAP.exitflag;
        end
        
        % Find non-dominated solutions
        nECs = problem.optimization.EC.num_per_dim;
        iECs = problem.optimization.EC.obj_num;
        fopt_min = min(fopt_anchor, [], 1);
        fopt_max = max(fopt_anchor, [], 1);
        % constraint combinations
        idxcon = ones(1, size(fopt_min, 2));
        idxcon(iECs) = 0;
        idxcon = logical(idxcon);
        epscon_arr = [];
        for idx = 1:size(fopt_min, 2)
            if idx ~= iECs
                tmp = transpose(linspace(fopt_min(idx), fopt_max(idx), nECs));
                epscon_arr{1, idx} = tmp(2:(end-1), :);
            else
                epscon_arr{1, idx} = fopt_max(1, idx);
            end
        end
        combinations = cell(1, numel(epscon_arr));
        [combinations{end:-1:1}] = ndgrid(epscon_arr{end:-1:1});
        combinations = cat(numel(epscon_arr)+1, combinations{:});
        combinations = reshape(combinations, [], numel(epscon_arr));
        ncomb = size(combinations, 1);
        % add random numbers to constraints
        for idx = 1:size(combinations, 2)
            if (idx ~= iECs) && problem.optimization.EC.randomize
                d = combinations(2,idx) - combinations(1,idx);
                combinations(:,idx) = combinations(:,idx) + d*(rand(size(combinations,1),1) - 0.5);
            end
        end
        % run optimization
        xopt_solution = zeros(ncomb, size(xopt_anchor, 2));
        fopt_solution = zeros(ncomb, size(fopt_anchor, 2));
        for idx = 1:ncomb
            epscon = combinations(idx, :);
            epscon_compact = repmat(epscon(1, idxcon), size(startptF, 1), 1);
            fprev_compact = startptF(:, idxcon);
            [~, idxcloser] = sortrows(sum((epscon_compact - fprev_compact).^2, 2), 1);
            xprevtmp = startptX(idxcloser, :);
            if size(xprevtmp, 1) == 0
                xprevtmp = (problem.bound.xlb + problem.bound.xub)./2;
            end
            ECprob.x0 = xprevtmp(1, :);
            if nMS < size(xprevtmp, 1)
                x0arr = xprevtmp(1:nMS, :);
            else
                x0arr = xprevtmp;
            end
            tpoints = CustomStartPointSet(x0arr);
            [xoptND, foptND] = runECs(problem, MS, ECprob, tpoints, ...
                iECs, epscon, surrF, surrC, surrCEQ, irmodel);
            xopt_solution(idx, :) = xoptND;
            fopt_solution(idx, :) = foptND;
        end
        xopt = [xopt_anchor; xopt_solution];
        fopt = [fopt_anchor; fopt_solution];
        out = [];
    otherwise
        error([solver, ' not supported']);
    end

    % Make solution unique
    num_x = size(xopt, 2);
    num_f = size(fopt, 2);
    xcomb = unique([xopt, fopt], 'rows');
    xP = xcomb(:, 1:num_x);
    fP = xcomb(:, (num_x + 1):(num_x + num_f));

    % Get constraints values
    [cP, ceqP] = surrogateCeval(xP, surrC, surrCEQ, irmodel);
    
    t_elapsed = toc;
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xopt, fopt, out] = runECs(problem, ECMS, ECprob, ECstartpts, ...
        nobj, epscon, surrF, surrC, surrCEQ, irmodel)
    xLast = [];
    myf = [];
    myc = [];
    myceq = [];
    fun = @ECsurrogateObj;
    cfun = @ECsurrogateNonlcon;
    
    ECprob.objective = @(x) ECsurrogateObj(x, surrF, surrC, surrCEQ, irmodel, nobj, epscon);
    ECprob.nonlcon = @(x) ECsurrogateNonlcon(x, surrF, surrC, surrCEQ, irmodel, nobj, epscon);

    [xopt, ~, exitflag, output, solutions] = run(ECMS, ECprob, ECstartpts);
    [fopt, ~] = surrogateFeval(xopt, surrF);
    out.exitflag = exitflag;
    out.output = output;
    out.solutions = solutions;
    
    if exitflag < 1
        warning(['MultiStart fmincon for epsilon-cpnstraints method does not converge.', ...
            ' Try problem.optimization.fmincon.solver = ''interior-point'' if not.']);
    end
    
    %----1---------2---------3---------4---------5---------6---------7---------8---------9---------0
    
    function f = ECsurrogateObj(x, surrF, surrC, surrCEQ, irmodel, nobj, epscon)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallECs(x, surrF, surrC, surrCEQ, irmodel);
            xLast = x;
        end
        f = myf(:, nobj);
    end
    
    %----1---------2---------3---------4---------5---------6---------7---------8---------9---------0
    
    function [c, ceq] = ECsurrogateNonlcon(x, surrF, surrC, surrCEQ, irmodel, nobj, epscon)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallECs(x, surrF, surrC, surrCEQ, irmodel);
            xLast = x;
        end
        if size(epscon, 1) == 0
            c = myc;
        else
            epsconmat = repmat(epscon, size(myf, 1), 1);
            ceps = myf - epsconmat;
            ceps(:, nobj) = 0;
            if size(myc, 1) == 0
                c = ceps;
            else
                c = [myc, ceps];
            end
        end
        ceq = myceq;
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [f, c, ceq] = computeallECs(x, surrF, surrC, surrCEQ, irmodel)
    [f, ~] = surrogateFeval(x, surrF);
    [c, ceq, ~, ~] = surrogateCeval(x, surrC, surrCEQ, irmodel);
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
