%% MO-ASMO-II :: compileStartingPoints function
% 1. If GA type solver, compile initial population from previous results
% 2. If epsilon-constraints type solver, choose initial points for the nonlinear program to start
% Usage:
%  ptout = compileStartingPoints(varargin)
% Arguments:
%  {problem, xhff, fhff, chff, ceqhff, xpred}
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function ptout = compileStartingPoints(problem, xhff, fhff, chff, ceqhff, xpred)
    solver = problem.optimization.solver;
    popsize = problem.optimization.nsga2.popsize;
    ineqtol = problem.sampling.tolerance.ineq;
    eqtol = problem.sampling.tolerance.eq;
    num_f = problem.bound.num_f;
    number = size(xhff, 1);

    c = max(chff, [], 2);
    ceq = max(ceqhff.^2, [], 2);
    
    if (size(c, 1) > 0) && (size(ceq, 1) > 0)
        idx_f = ((c - ineqtol) <= 0) & ((ceq - eqtol) <= 0);
    elseif (size(c, 1) > 0)
        idx_f = ((c - ineqtol) <= 0);
    elseif (size(ceq, 1) > 0)
        idx_f = ((ceq - eqtol) <= 0);
    else
        idx_f = [];
    end
    xhff_f = xhff(idx_f, :);
    fhff_f = fhff(idx_f, :);
    %chff_f = chff(idx_f); % >> to check if c is satisfied. to be deleted later.
    %ceqhff_f = ceqhff(idx_f); % >> to check if ceq is satisfied. to be deleted later.
    number_f = size(xhff_f, 1);

    ptout = zeros(0, size(xhff,2));

    switch lower(solver)
    case 'nsga-ii'
        if (size(xhff, 1) < 20) % if high fidelity result pool is too small
            ptout = xhff;
        else
            % Sample starting points from high fidelity results, which satisfies constraints
            if number_f > 0
                if number_f < 20
                    ptout = xhff_f;
                else
                    [xhff_f_sort, ~, ihff_f_sort] = ndSort(xhff_f, fhff_f);
                    if (max(ihff_f_sort) == 1) % has only one rank
                        pttmp = xhff_f_sort((ihff_f_sort == 1), :);
                    else
                        for idx = max(ihff_f_sort):(-1):1
                            pttmp = xhff_f_sort((ihff_f_sort <= idx), :);
                            if (size(pttmp, 1) < 0.3*popsize), break; end
                        end
                    end
                    ptout = vertcat(ptout, pttmp);
                end
            end
            % Sample starting points from high fidelity results, regardless constraint compliance
            if number < 20
                ptout = xhff;
            else
                [xhff_sort, ~, ihff_sort] = ndSort(xhff, fhff);
                if (max(ihff_sort) == 1) % has only one rank
                    pttmp = xhff_sort((ihff_sort == 1), :);
                else
                    for idx = max(ihff_sort):(-1):1
                        pttmp = xhff_sort((ihff_sort <= idx), :);
                        if (size(pttmp, 1) < 0.3*popsize), break; end
                    end
                end
                ptout = vertcat(ptout, pttmp);
            end
            ptout = unique(ptout, 'rows');
            % Sampling starting points from previously predicted points
            nreq = size(ptout, 1);
            if nreq < popsize
                nreq = popsize - nreq;
                if nreq > size(xpred, 1)
                    ptout = vertcat(ptout, xpred);
                else
                    pttmp = datasample(xpred, nreq, 'Replace', false);
                    ptout = vertcat(ptout, pttmp);
                end
            end
        end
    case 'epsilon-constraints'
        [xhff_f_sort, ~, ihff_f_sort] = ndSort(xhff_f, fhff_f);
        xhff_f_nd = xhff_f_sort(ihff_f_sort == 1);
        ptout = vertcat(ptout, xhff_f_nd);
    otherwise
        error([solver, ' not supported']);
    end

    % Choose unique points
    ptout = unique(ptout, 'rows');
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
