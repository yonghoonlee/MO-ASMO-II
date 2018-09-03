%% MO-ASMO-II :: compileStartingPoints function
% 1. If GA type solver, compile initial population from previous results
% 2. If epsilon-constraints type solver, choose initial points for the nonlinear program to start
% Usage:
%  pts = compileStartingPoints(varargin)
% Arguments:
%  {problem, xhff, fhff, chff, ceqhff, xpred, fpred}
% Returns:
%  pts.{x, f}
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function pts = compileStartingPoints(problem, xhff, fhff, chff, ceqhff, xpred, fpred)
    solver = problem.optimization.solver;
    popsize = problem.optimization.nsga2.popsize;
    ineqtol = problem.sampling.tolerance.ineq;
    eqtol = problem.sampling.tolerance.eq;
    num_f = problem.bound.num_f;
    number = size(xhff, 1);

    c = max(chff, [], 2);
    ceq = max(ceqhff.^2, [], 2);
    
    if (size(c, 2) > 0) && (size(ceq, 2) > 0)
        idx_f = ((c - ineqtol) <= 0) & ((ceq - eqtol) <= 0);
    elseif (size(c, 2) > 0)
        idx_f = ((c - ineqtol) <= 0);
    elseif (size(ceq, 2) > 0)
        idx_f = ((ceq - eqtol) <= 0);
    else
        idx_f = [];
    end
    xhff_f = xhff(idx_f, :);
    fhff_f = fhff(idx_f, :);
    number_f = size(xhff_f, 1);

    ptoutX = zeros(0, size(xhff,2));
    ptoutF = zeros(0, size(fhff,2));

    switch lower(solver)
    case 'nsga-ii'
        if (number < 20) % if high fidelity result pool is too small
            ptoutX = xhff;
            ptoutF = fhff;
        else
            % Sample starting points from high fidelity results, which satisfies constraints
            if number_f > 0
                if number_f < 20
                    ptoutX = xhff_f;
                    ptoutF = fhff_f;
                else
                    [xhff_f_sort, fhff_f_sort, ihff_f_sort] = ndSort(xhff_f, fhff_f);
                    if (max(ihff_f_sort) == 1) % has only one rank
                        pttmpX = xhff_f_sort((ihff_f_sort == 1), :);
                        pttmpF = fhff_f_sort((ihff_f_sort == 1), :);
                    else
                        for idx = max(ihff_f_sort):(-1):1
                            pttmpX = xhff_f_sort((ihff_f_sort <= idx), :);
                            pttmpF = fhff_f_sort((ihff_f_sort <= idx), :);
                            if (size(pttmpX, 1) < 0.3*popsize), break; end
                        end
                    end
                    ptoutX = vertcat(ptoutX, pttmpX);
                    ptoutF = vertcat(ptoutF, pttmpF);
                end
            end
            % Sample starting points from high fidelity results, regardless constraint compliance
            [xhff_sort, fhff_sort, ihff_sort] = ndSort(xhff, fhff);
            if (max(ihff_sort) == 1) % has only one rank
                pttmpX = xhff_sort((ihff_sort == 1), :);
                pttmpF = fhff_sort((ihff_sort == 1), :);
            else
                for idx = max(ihff_sort):(-1):1
                    pttmpX = xhff_sort((ihff_sort <= idx), :);
                    pttmpF = fhff_sort((ihff_sort <= idx), :);
                    if (size(pttmpX, 1) < 0.3*popsize), break; end
                end
            end
            ptoutX = vertcat(ptoutX, pttmpX);
            ptoutF = vertcat(ptoutF, pttmpF);
            [~, ia, ~] = unique(ptoutX, 'rows');
            ptoutX = ptoutX(ia, :);
            ptoutF = ptoutF(ia, :);
            % Sampling starting points from previously predicted points
            nreq = size(ptoutX, 1);
            if nreq < popsize
                nreq = popsize - nreq;
                if nreq > size(xpred, 1)
                    ptoutX = vertcat(ptoutX, xpred);
                    ptoutF = vertcat(ptoutF, fpred);
                else
                    [~, ia] = datasample(xpred, nreq, 'Replace', false);
                    pttmpX = xpred(ia, :);
                    pttmpF = fpred(ia, :);
                    ptoutX = vertcat(ptoutX, pttmpX);
                    ptoutF = vertcat(ptoutF, pttmpF);
                end
            end
        end
    case 'epsilon-constraints'
        [xhff_f_sort, fhff_f_sort, ihff_f_sort] = ndSort(xhff_f, fhff_f);
        xhff_f_nd = xhff_f_sort(ihff_f_sort == 1, :);
        fhff_f_nd = fhff_f_sort(ihff_f_sort == 1, :);
        ptoutX = vertcat(ptoutX, xhff_f_nd);
        ptoutF = vertcat(ptoutF, fhff_f_nd);
    otherwise
        error([solver, ' not supported']);
    end

    % Choose unique points
    [~, ia, ~] = unique(ptoutX, 'rows');
    ptoutX = ptoutX(ia, :);
    ptoutF = ptoutF(ia, :);
    
    % Return
    pts.x = ptoutX;
    pts.f = ptoutF;
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
