%% MO-ASMO-II :: debugAnalysis function
% 1. Plot and display debug information
% Usage:
%  debugAnalysis(problem, subroutine_name, varargin)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function debugAnalysis(problem, subroutine_name, varargin)
    declareGlobalVariables;
    declareDebugVariables;
    
    if verbose == 2
        switch lower(subroutine_name)
        case 'initialsampling'
            if (nargin < 3), error('Debugging initial sampling requires xt'); end
            xt = varargin{1};
            xt = varScale(xt, problem.bound.xlb, problem.bound.xub, 'scale');
            disp(['Number of initial samples created: ', num2str(size(xt, 1))]);
            if size(xt, 1) == 0
                warning('No initial samples generated');
                return;
            end
            
            % Design variable space
            if size(xt, 1) == 0, return; end
            fg_debug1 = debug_figure_open(fg_debug1);
            phA = []; lhA = {};
            if (size(xt, 2) == 1)
                xtx = 1:size(xt, 1);
                hold off;
                ph1 = plot(xtx, reshape(xt, 1, numel(xt)), 'rs'); hold on;
                phA = horzcat(phA, ph1);
                lhA = horzcat(lhA, {'feasible point in $x$-space'});
                xlabel('Training sample number');
                ylabel('Normalized $x$');
            elseif (size(xt, 2) == 2)
                hold off;
                ph1 = plot(xt(:,1), xt(:,2), 'rs'); hold on;
                phA = horzcat(phA, ph1);
                lhA = horzcat(lhA, {'feasible point in $x$-space'});
                xlabel('Normalized $x_1$');
                xlabel('Normalized $x_2$');
            else
                xtx = 1:size(xt, 2);
                cm = plasma(size(xt, 1));
                hold off;
                for idx = 1:size(xt, 1)
                    xty = xt(idx, :);
                    ph1 = plot(xtx, xty, '-', 'Color', cm(idx, :)); hold on;
                end
                phA = horzcat(phA, ph1);
                lhA = horzcat(lhA, {'feasible point in $x$-space'});
                xlabel('Variable number ($\#$)');
                ylabel('Normalized value ($x_\#$)');
            end
            legend(phA, lhA, 'Location', 'Best');
            
            % Objective function space
            if problem.functions.hifi_expensive, return; end
            xt = varScale(xt, problem.bound.xlb, problem.bound.xub, 'unscale');
            [hffF, hffC, hffCEQ] = evaluateHff(problem, xt);
            [idxValid, ~] = separateNaN(xt, hffF, hffC, hffCEQ);
            xt = xt(idxValid, :); hffF = hffF(idxValid, :);
            if (size(hffC, 2) ~= 0), hffC = hffC(idxValid, :);
            else hffC = zeros(size(hffF, 1), 1); end
            hffC = max(hffC, [], 2);
            if (size(hffCEQ, 2) ~= 0), hffCEQ = hffCEQ(idxValid, :);
            else hffCEQ = zeros(size(hffF, 1), 1); end
            hffCEQ = sqrt(max(hffCEQ.^2, [], 2));
            ioptm = (((hffC - problem.control.tolC) <= 0) ...
                & ((hffCEQ - problem.control.tolCEQ) <= 0));
            ioptm = enforceIndexLincon(problem, ioptm, xt);
            xts = xt(ioptm, :); xfs = hffF(ioptm, :);
            xtn = xt(~ioptm, :); xfn = hffF(~ioptm, :);
            [~, xfs, xis] = ndSort(xts, xfs);
            [~, xfn, ~] = ndSort(xtn, xfn);
            fg_debug2 = debug_figure_open(fg_debug2);
            phB = []; lhB = {};
            % Feasible points
            if (size(xfs, 1) ~= 0)
                if (size(xfs, 2) == 1)
                    xfx = 1:size(xfs, 1);
                    hold off;
                    ph2 = plot(xfx, reshape(xfs, 1, numel(xfs)), 'rs');
                    phB = horzcat(phB, ph2);
                    lhB = horzcat(lhB, {'feasible point in $f$-space'});
                    xlabel('Training sample number');
                    ylabel('Normalized $f$');
                elseif (size(xfs, 2) == 2)
                    hold off;
                    ph2 = plot(xfs(:,1), xfs(:,2), 'rs'); hold on;
                    phB = horzcat(phB, ph2);
                    lhB = horzcat(lhB, {'feasible point in $f$-space'});
                    xfss = sortrows(xfs(xis == 1, :), 1);
                    ph3 = plot(xfss(:, 1), xfss(:, 2), 'r-'); hold on;
                    phB = horzcat(phB, ph3);
                    lhB = horzcat(lhB, {'non-dominated point'});
                    xlabel('Normalized $f_1$');
                    xlabel('Normalized $f_2$');
                else
                    xfx = 1:size(xfs, 2);
                    cm = plasma(max(size(xfs, 1), 1));
                    hold off;
                    for idx = 1:size(xfs, 1)
                        xfy = xfs(idx, :);
                        ph2 = plot(xfx, xfy, '-', 'Color', cm(idx, :)); hold on;
                    end
                    phB = horzcat(phB, ph2);
                    lhB = horzcat(lhB, {'feasible point in $f$-space'});
                    xlabel('Variable number ($\#$)');
                    ylabel('Normalized value ($f_\#$)');
                end
            end
            % Infeasible points
            if (size(xfn, 1) ~= 0)
                if (size(xfn, 2) == 1)
                    xfx = 1:size(xfn, 1);
                    ph4 = plot(xfx, reshape(xfn, 1, numel(xfn)), 'ms');
                    phB = horzcat(phB, ph4);
                    lhB = horzcat(lhB, {'infeasible point in $f$-space'});
                elseif (size(xfn, 2) == 2)
                    ph4 = plot(xfn(:,1), xfn(:,2), 'ms'); hold on;
                    phB = horzcat(phB, ph4);
                    lhB = horzcat(lhB, {'infeasible point in $f$-space'});
                else
                    xfx = 1:size(xfn, 2);
                    for idx = 1:size(xfn, 1)
                        xfy = xfn(idx, :);
                        ph4 = plot(xfx, xfy, 'm-'); hold on;
                    end
                    phB = horzcat(phB, ph4);
                    lhB = horzcat(lhB, {'infeasible point in $f$-space'});
                end
            end
            legend(phB, lhB, 'Location', 'Best');
            
        case 'updatesampling'
            if (nargin < 5), error('Debugging update sampling requires xt1, xt2, xP'); end
            xt1 = varargin{1};
            xt2 = varargin{2};
            xP = varargin{3};
            xt1 = varScale(xt1, problem.bound.xlb, problem.bound.xub, 'scale');
            xt2 = varScale(xt2, problem.bound.xlb, problem.bound.xub, 'scale');
            xP = varScale(xP, problem.bound.xlb, problem.bound.xub, 'scale');
            disp(['Number of update samples created: ', num2str(size([xt1; xt2], 1))]);
            disp(['  Number of exploration samples:  ', num2str(size(xt1, 1))]);
            disp(['  Number of exploitation samples: ', num2str(size(xt2, 1))]);
            if size([xt1; xt2], 1) == 0
                warning('No update samples generated');
                return;
            end
            
            % Design variable space
            if ((size(xt1, 1) == 0) && (size(xt2, 1) == 0)), return; end
            fg_debug1 = debug_figure_open(fg_debug1);
            phC = []; lhC = {};
            % Pareto set
            if (size(xP, 1) ~= 0)
                if (size(xP, 2) == 1)
                    xtx = 1:size(xP, 1);
                    hold off;
                    ph1 = plot(xtx, reshape(xP, 1, numel(xP)), 'gh'); hold on;
                    phC = horzcat(phC, ph1);
                    lhC = horzcat(lhC, {'predicted Pareto set'});
                    xlabel('Training sample number');
                    ylabel('Normalized $x$');
                elseif (size(xP, 2) == 2)
                    hold off;
                    ph1 = plot(xP(:,1), xP(:,2), 'gh'); hold on;
                    phC = horzcat(phC, ph1);
                    lhC = horzcat(lhC, {'predicted Pareto set'});
                    xlabel('Normalized $x_1$');
                    xlabel('Normalized $x_2$');
                else
                    xtx = 1:size(xP, 2);
                    hold off;
                    for idx = 1:size(xP, 1)
                        xty = xP(idx, :);
                        ph1 = plot(xtx, xty, 'g-'); hold on;
                    end
                    phC = horzcat(phC, ph1);
                    lhC = horzcat(lhC, {'predicted Pareto set'});
                    xlabel('Variable number ($\#$)');
                    ylabel('Normalized value ($x_\#$)');
                end
            end
            % Feasible exploration samples
            if (size(xt1, 1) ~= 0)
                if (size(xt1, 2) == 1)
                    xtx = 1:size(xt1, 1);
                    ph2 = plot(xtx, reshape(xt1, 1, numel(xt1)), 'rs'); hold on;
                    phC = horzcat(phC, ph2);
                    lhC = horzcat(lhC, {'exploration point in $x$-space'});
                elseif (size(xt1, 2) == 2)
                    ph2 = plot(xt1(:,1), xt1(:,2), 'rs'); hold on;
                    phC = horzcat(phC, ph2);
                    lhC = horzcat(lhC, {'exploration point in $x$-space'});
                else
                    xtx = 1:size(xt1, 2);
                    cm = plasma(ceil(1.2*(size(xt1, 1) + size(xt2, 1))));
                    for idx = 1:size(xt1, 1)
                        xty1 = xt1(idx, :);
                        ph2 = plot(xtx, xty1, '-', 'Color', cm(idx, :)); hold on;
                    end
                    phC = horzcat(phC, ph2);
                    lhC = horzcat(lhC, {'exploration point in $x$-space'});
                end
            end
            % Feasible exploitation samples
            if (size(xt2, 1) ~= 0)
                if (size(xt2, 2) == 1)
                    xtx = 1:size(xt2, 1);
                    ph3 = plot(xtx, reshape(xt2, 1, numel(xt2)), 'bd'); hold on;
                    phC = horzcat(phC, ph3);
                    lhC = horzcat(lhC, {'exploitation point in $x$-space'});
                elseif (size(xt2, 2) == 2)
                    ph3 = plot(xt2(:,1), xt2(:,2), 'bd'); hold on;
                    phC = horzcat(phC, ph3);
                    lhC = horzcat(lhC, {'exploitation point in $x$-space'});
                else
                    xtx = 1:size(xt2, 2);
                    for idx = 1:size(xt2, 1)
                        xty2 = xt2(idx, :);
                        ph3 = plot(xtx, xty2, '-', 'Color', cm(end-idx+1, :)); hold on;
                    end
                    phC = horzcat(phC, ph3);
                    lhC = horzcat(lhC, {'exploitation point in $x$-space'});
                end
            end
            legend(phC, lhC, 'Location', 'Best');
            
            % Objective function space
            if problem.functions.hifi_expensive, return; end
            % Exploration samples
            xt = xt1;
            xt = varScale(xt, problem.bound.xlb, problem.bound.xub, 'unscale');
            [hffF, hffC, hffCEQ] = evaluateHff(problem, xt);
            [idxValid, ~] = separateNaN(xt, hffF, hffC, hffCEQ);
            xt = xt(idxValid, :); hffF = hffF(idxValid, :);
            if (size(hffC, 2) ~= 0), hffC = hffC(idxValid, :);
            else hffC = zeros(size(hffF, 1), 1); end
            hffC = max(hffC, [], 2);
            if (size(hffCEQ, 2) ~= 0), hffCEQ = hffCEQ(idxValid, :);
            else hffCEQ = zeros(size(hffF, 1), 1); end
            hffCEQ = sqrt(max(hffCEQ.^2, [], 2));
            ioptm = (((hffC - problem.control.tolC) <= 0) ...
                & ((hffCEQ - problem.control.tolCEQ) <= 0));
            ioptm = enforceIndexLincon(problem, ioptm, xt);
            xts = xt(ioptm, :); xfs = hffF(ioptm, :);
            xtn = xt(~ioptm, :); xfn = hffF(~ioptm, :);
            [~, xfs, xis] = ndSort(xts, xfs);
            [~, xfn, ~] = ndSort(xtn, xfn);
            fg_debug2 = debug_figure_open(fg_debug2);
            phD = []; lhD = {};
            % Feasible exploration samples
            if (size(xfs, 1) ~= 0)
                if (size(xfs, 2) == 1)
                    xfx = 1:size(xfs, 1);
                    hold off;
                    ph4 = plot(xfx, reshape(xfs, 1, numel(xfs)), 'rs');
                    phD = horzcat(phD, ph4);
                    lhD = horzcat(lhD, {'feasible exploration point in $f$-space'});
                    xlabel('Training sample number');
                    ylabel('Normalized $f$');
                elseif (size(xfs, 2) == 2)
                    hold off;
                    ph4 = plot(xfs(:,1), xfs(:,2), 'rs'); hold on;
                    phD = horzcat(phD, ph4);
                    lhD = horzcat(lhD, {'feasible exploration point in $f$-space'});
                    xfss = sortrows(xfs(xis == 1, :), 1);
                    ph5 = plot(xfss(:, 1), xfss(:, 2), 'r-'); hold on;
                    phD = horzcat(phD, ph5);
                    lhD = horzcat(lhD, {'non-dominated point'});
                    xlabel('Normalized $f_1$');
                    xlabel('Normalized $f_2$');
                else
                    xfx = 1:size(xfs, 2);
                    cm = plasma(max(ceil(1.2*size([xt1;xt2], 1)), 1));
                    hold off;
                    for idx = 1:size(xfs, 1)
                        xfy = xfs(idx, :);
                        ph4 = plot(xfx, xfy, '-', 'Color', cm(idx, :)); hold on;
                    end
                    phD = horzcat(phD, ph4);
                    lhD = horzcat(lhD, {'feasible exploration point in $f$-space'});
                    xlabel('Variable number ($\#$)');
                    ylabel('Normalized value ($f_\#$)');
                end
            end
            % Ineasible exploration samples
            if (size(xfn, 1) ~= 0)
                if (size(xfn, 2) == 1)
                    xfx = 1:size(xfn, 1)
                    ph6 = plot(xfx, reshape(xfn, 1, numel(xfn)), 'ms');
                    phD = horzcat(phD, ph6);
                    lhD = horzcat(lhD, {'infeasible exploration point in $f$-space'});
                elseif (size(xfn, 2) == 2)
                    ph6 = plot(xfn(:,1), xfn(:,2), 'ms'); hold on;
                    phD = horzcat(phD, ph6);
                    lhD = horzcat(lhD, {'infeasible exploration point in $f$-space'});
                else
                    xfx = 1:size(xfn, 2);
                    for idx = 1:size(xfn, 1)
                        xfy = xfn(idx, :);
                        ph6 = plot(xfx, xfy, 'm-'); hold on;
                    end
                    phD = horzcat(phD, ph6);
                    lhD = horzcat(lhD, {'infeasible exploration point in $f$-space'});
                end
            end
            % Exploitation samples
            xt = xt2;
            xt = varScale(xt, problem.bound.xlb, problem.bound.xub, 'unscale');
            [hffF, hffC, hffCEQ] = evaluateHff(problem, xt);
            [idxValid, ~] = separateNaN(xt, hffF, hffC, hffCEQ);
            xt = xt(idxValid, :); hffF = hffF(idxValid, :);
            if (size(hffC, 2) ~= 0), hffC = hffC(idxValid, :);
            else hffC = zeros(size(hffF, 1), 1); end
            hffC = max(hffC, [], 2);
            if (size(hffCEQ, 2) ~= 0), hffCEQ = hffCEQ(idxValid, :);
            else hffCEQ = zeros(size(hffF, 1), 1); end
            hffCEQ = sqrt(max(hffCEQ.^2, [], 2));
            ioptm = (((hffC - problem.control.tolC) <= 0) ...
                & ((hffCEQ - problem.control.tolCEQ) <= 0));
            ioptm = enforceIndexLincon(problem, ioptm, xt);
            xts = xt(ioptm, :); xfs = hffF(ioptm, :);
            xtn = xt(~ioptm, :); xfn = hffF(~ioptm, :);
            [~, xfs, xis] = ndSort(xts, xfs);
            [~, xfn, xin] = ndSort(xtn, xfn);
            fg_debug2 = debug_figure_open(fg_debug2);
            % Feasible exploitation samples
            if (size(xfs, 1) ~= 0)
                if (size(xfs, 2) == 1)
                    xfx = 1:size(xfs, 1);
                    ph7 = plot(xfx, reshape(xfs, 1, numel(xfs)), 'bd');
                    phD = horzcat(phD, ph7);
                    lhD = horzcat(lhD, {'feasible exploitation point in $f$-space'});
                elseif (size(xfs, 2) == 2)
                    ph7 = plot(xfs(:,1), xfs(:,2), 'bd'); hold on;
                    phD = horzcat(phD, ph7);
                    lhD = horzcat(lhD, {'feasible exploitation point in $f$-space'});
                    xfss = sortrows(xfs(xis == 1, :), 1);
                    ph8 = plot(xfss(:, 1), xfss(:, 2), 'b-'); hold on;
                    phD = horzcat(phD, ph8);
                    lhD = horzcat(lhD, {'non-dominated point'});
                else
                    xfx = 1:size(xfs, 2);
                    for idx = 1:size(xfs, 1)
                        xfy = xfs(idx, :);
                        ph7 = plot(xfx, xfy, '-', 'Color', cm(end-idx+1, :)); hold on;
                    end
                    phD = horzcat(phD, ph7);
                    lhD = horzcat(lhD, {'feasible exploitation point in $f$-space'});
                end
            end
            % Infeasible exploitation samples
            if (size(xfn, 1) ~= 0)
                if (size(xfn, 2) == 1)
                    xfx = 1:size(xfn, 1);
                    ph9 = plot(xfx, reshape(xfn, 1, numel(xfn)), 'cd');
                    phD = horzcat(phD, ph9);
                    lhD = horzcat(lhD, {'infeasible exploitation point in $f$-space'});
                elseif (size(xfn, 2) == 2)
                    ph9 = plot(xfn(:,1), xfn(:,2), 'cd'); hold on;
                    phD = horzcat(phD, ph9);
                    lhD = horzcat(lhD, {'infeasible exploitation point in $f$-space'});
                else
                    xfx = 1:size(xfn, 2);
                    for idx = 1:size(xfn, 1)
                        xfy = xfn(idx, :);
                        ph9 = plot(xfx, xfy, 'c-'); hold on;
                    end
                    phD = horzcat(phD, ph9);
                    lhD = horzcat(lhD, {'infeasible exploitation point in $f$-space'});
                end
            end
            % Pareto set
            xt = xP;
            xt = varScale(xt, problem.bound.xlb, problem.bound.xub, 'unscale');
            [hffF, hffC, hffCEQ] = evaluateHff(problem, xt);
            [idxValid, ~] = separateNaN(xt, hffF, hffC, hffCEQ);
            xt = xt(idxValid, :); hffF = hffF(idxValid, :);
            if (size(hffC, 2) ~= 0), hffC = hffC(idxValid, :);
            else hffC = zeros(size(hffF, 1), 1); end
            hffC = max(hffC, [], 2);
            if (size(hffCEQ, 2) ~= 0), hffCEQ = hffCEQ(idxValid, :);
            else hffCEQ = zeros(size(hffF, 1), 1); end
            hffCEQ = sqrt(max(hffCEQ.^2, [], 2));
            ioptm = (((hffC - problem.control.tolC) <= 0) ...
                & ((hffCEQ - problem.control.tolCEQ) <= 0));
            ioptm = enforceIndexLincon(problem, ioptm, xt);
            xfs = hffF(ioptm, :);
            xfs = sortrows(xfs, 1);
            fg_debug2 = debug_figure_open(fg_debug2);
            if (size(xfs, 1) ~= 0)
                if (size(xfs, 2) == 1)
                    xfx = 1:size(xfs, 1)
                    ph10 = plot(xfx, reshape(xfs, 1, numel(xfs)), 'gh');
                    phD = horzcat(phD, ph10);
                    lhD = horzcat(lhD, {'predicted Pareto set in $f$-space'});
                elseif (size(xfs, 2) == 2)
                    ph10 = plot(xfs(:,1), xfs(:,2), 'gh'); hold on;
                    phD = horzcat(phD, ph10);
                    lhD = horzcat(lhD, {'predicted Pareto set in $f$-space'});
                else
                    xfx = 1:size(xfs, 2);
                    for idx = 1:size(xfs, 1)
                        xfy = xfs(idx, :);
                        ph10 = plot(xfx, xfy, 'g-'); hold on;
                    end
                    phD = horzcat(phD, ph10);
                    lhD = horzcat(lhD, {'predicted Pareto set in $f$-space'});
                end
            end
            legend(phD, lhD, 'Location', 'Best');
            
        case 'force-directed'
            if (nargin < 5), error('Debugging force-directed sampling requires xM, xB, xP'); end
            xM = varargin{1};
            xB = varargin{2};
            xP = varargin{3};
            if (size(xM, 1) == 0) || (size(xB, 1) == 0) || (size(xP, 1) == 0), return; end
            
            fg_debug3 = debug_figure_open(fg_debug3);
            % Pareto set
            if (size(xP, 1) ~= 0)
                if (size(xP, 2) == 1)
                    xtx = 1:size(xP, 1);
                    xty = reshape(xP, 1, numel(xP));
                    hold off;
                    plot(xtx, xty, 'g.'); hold on;
                elseif (size(xP, 2) == 2)
                    hold off;
                    plot(xP(:, 1), xP(:, 2), 'g.'); hold on;
                elseif (size(xP, 2) == 3)
                    hold off;
                    plot3(xP(:, 1), xP(:, 2), xP(:, 3), 'g.'); hold on;
                elseif (size(xP, 2) == 6)
                    subplot(1,2,1); hold off;
                    plot3(xP(:, 1), xP(:, 2), xP(:, 3), 'g.'); hold on;
                    subplot(1,2,2); hold off;
                    plot3(xP(:, 4), xP(:, 5), xP(:, 6), 'g.'); hold on;
                else
                    hold off;
                    xtx = 1:size(xP, 2);
                    for idx = 1:size(xP, 1)
                        plot(xtx, xP(idx, :), 'g-'); hold on;
                    end
                end
            end
            % Base points
            if (size(xB, 1) ~= 0)
                if (size(xB, 2) == 1)
                    xtx = 1:size(xB, 1);
                    xty = reshape(xB, 1, numel(xB));
                    plot(xtx, xty, 'bo'); hold on;
                elseif (size(xB, 2) == 2)
                    plot(xB(:, 1), xB(:, 2), 'bo'); hold on;
                elseif (size(xB, 2) == 3)
                    plot3(xB(:, 1), xB(:, 2), xB(:, 3), 'bo'); hold on;
                elseif (size(xB, 2) == 6)
                    subplot(1,2,1);
                    plot3(xB(:, 1), xB(:, 2), xB(:, 3), 'bo'); hold on;
                    subplot(1,2,2);
                    plot3(xB(:, 4), xB(:, 5), xB(:, 6), 'bo'); hold on;
                else
                    xtx = 1:size(xB, 2);
                    for idx = 1:size(xB, 1)
                        plot(xtx, xB(idx, :), 'b-'); hold on;
                    end
                end
            end
            % Moving points
            if (size(xM, 1) ~= 0)
                if (size(xM, 2) == 1)
                    xtx = 1:size(xM, 1);
                    xty = reshape(xM, 1, numel(xM));
                    plot(xtx, xty, 'rh'); hold on;
                elseif (size(xM, 2) == 2)
                    plot(xM(:, 1), xM(:, 2), 'rh'); hold on;
                elseif (size(xM, 2) == 3)
                    plot3(xM(:, 1), xM(:, 2), xM(:, 3), 'rh'); hold on;
                elseif (size(xM, 2) == 6)
                    subplot(1,2,1);
                    plot3(xM(:, 1), xM(:, 2), xM(:, 3), 'rh'); hold on;
                    subplot(1,2,2);
                    plot3(xM(:, 4), xM(:, 5), xM(:, 6), 'rh'); hold on;
                else
                    xtx = 1:size(xM, 2);
                    for idx = 1:size(xM, 1)
                        plot(xtx, xM(idx, :), 'r-'); hold on;
                    end
                end
            end
            ax = gca; limits = axis(ax);
            pause(0.5);
            
        otherwise
        end
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function fg = debug_figure_open(fg)
    try
        figure(fg);
        fg.Visible = 'on';
    catch
        fg = figure('Color',[1 1 1]);
        fg.Visible = 'on';
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
