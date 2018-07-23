%% MO-ASMO-II :: surrogateOptim function
% 1. Using surrogate models for f, c, ceq, run surrogate-based optimization
% Usage:
%  varargout = surrogateOptim(varargin)
% Arguments:
%  {problem, surrF, surrC, surrCEQ, irmodel, startpts}
% Returns:
%  {xP, fP, cP, ceqP, out}
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xP, fP, cP, ceqP, out] = surrogateOptim(problem, surrF, surrC, surrCEQ, irmodel, startpts)
    declareGlobalVariables;
    solver = problem.optimization.solver;
    if verbose, disp('Solving multiobjective optimization using surrogate models'); end

    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    lb = reshape(problem.bound.xlb, 1, problem.bound.num_x);
    ub = reshape(problem.bound.xub, 1, problem.bound.num_x);

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

        opt = gaoptimset(@gamultiobj);
        opt.Vectorized = 'on';
        opt.UseParallel = false;
        opt.PopulationSize = popsize;
        opt.ParetoFraction = paretofrac;
        opt.StallGenLimit = stallgenlimit;
        opt.InitialPopulation = startpts;
        opt.PlotFcns = @gaplotpareto;

        [xopt, fopt, exitflag, output, population, score] = gamultiobj( ...
            @(x) surrogateFeval(x, surrF), nxvar, A, b, Aeq, beq, lb, ub, ...
            @(x) surrogateCeval(x, surrC, surrCEQ), opt);
        
        out.exitflag = exitflag;
        out.output = output;
    case 'epsilon-constraints'






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
    [cP, ceqP] = surrogateCeval(xP, surrC, surrCEQ);
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
