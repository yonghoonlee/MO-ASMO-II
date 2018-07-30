%% MO-ASMO-II :: debugAnalysis function
% 1. Plot and display debug information
% Usage:
%  debugAnalysis(subroutine, varargin)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function debugAnalysis(problem, subroutine, varargin)
    declareGlobalVariables;
    declareDebugVariables;
    
    if verbose == 2
        switch lower(subroutine)
        case 'initialsampling'
            if (nargin < 3), error('Debugging initial sampling requires xt'); end
            xt = varargin{1};
            xt = varScale(xt, problem.bound.xlb, problem.bound.xub, 'scale');
            disp(['Number of initial samples created: ', num2str(size(xt, 1))]);
            % Design variable space
            if size(xt, 1) == 0, return; end
            fg_debug1 = debug_figure_open(fg_debug1);
            if (size(xt, 2) == 1)
                xtx = 1:size(xt, 1);
                hold off;
                plot(xtx, reshape(xt, 1, numel(xt)), 'rs'); hold on;
                xlabel('Training sample number');
                ylabel('Normalized $x$');
            elseif (size(xt, 2) == 2)
                hold off;
                plot(xt(:,1), xt(:,2), 'rs'); hold on;
                xlabel('Normalized $x_1$');
                xlabel('Normalized $x_2$');
            else
                xtx = 1:size(xt, 2);
                cm = plasma(size(xt, 1));
                hold off;
                for idx = 1:size(xt, 1)
                    xty = xt(idx, :);
                    plot(xtx, xty, '-', 'Color', cm(idx, :)); hold on;
                end
                xlabel('Variable number ($\#$)');
                ylabel('Normalized value ($x_\#$)');
            end
            % Objective function space
            if problem.functions.hifi_expensive, return; end
            [hffF, hffC, hffCEQ] = evaluateHff(problem, xt);
            [idxValid, ~] = separateNaN(xt, hffF, hffC, hffCEQ);
            xt = xt(idxValid, :); hffF = hffF(idxValid, :);
            if (size(hffC, 1) ~= 0), hffC = hffC(idxValid, :);
            else hffC = zeros(size(hffF, 1), 1); end
            hffC = max(hffC, [], 2);
            if (size(hffCEQ, 1) ~= 0), hffCEQ = hffCEQ(idxValid, :);
            else hffCEQ = zeros(size(hffCEQ, 1), 1); end
            hffCEQ = sqrt(max(hffCEQ.^2, [], 2));
            ioptm = (((hffC - problem.control.tolC) <= 0) ...
                & ((hffCEQ - problem.control.tolCEQ) <= 0));
            xts = xt(ioptm, :); xfs = hffF(ioptm, :);
            [xts, xfs, xis] = ndSort(xts, xfs);
            fg_debug2 = debug_figure_open(fg_debug2);
            if (size(xfs, 2) == 1)
                xfx = 1:size(xfs, 1)
                hold off;
                plot(xfx, reshape(xfs, 1, numel(xfs)), 'rs');
                xlabel('Training sample number');
                ylabel('Normalized $f$');
            elseif (size(xfs, 2) == 2)
                hold off;
                plot(xfs(:,1), xfs(:,2), 'rs'); hold on;
                plot(xfs(xis == 1, 1), xfs(xis == 1, 2), 'r-'); hold on;
                xlabel('Normalized $f_1$');
                xlabel('Normalized $f_2$');
            else
                xfx = 1:size(xfs, 2);
                cm = plasma(max(size(xfs, 1), 1));
                hold off;
                for idx = 1:size(xfs, 1)
                    xfy = xfs(idx, :);
                    plot(xfx, xfy, '-', 'Color', cm(idx, :)); hold on;
                end
                xlabel('Variable number ($\#$)');
                ylabel('Normalized value ($f_\#$)');
            end
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
            % Design variable space
            if ((size(xt1, 1) == 0) && (size(xt2, 1) == 0)), return; end
            fg_debug1 = debug_figure_open(fg_debug1);
            if (size(xt1, 2) == 1)
                xtx = 1:size(xt1, 1);
                hold off;
                plot(xtx, reshape(xt1, 1, numel(xt1)), 'rs'); hold on;
                xlabel('Training sample number');
                ylabel('Normalized $x$');
            elseif (size(xt1, 2) == 2)
                hold off;
                plot(xt1(:,1), xt1(:,2), 'rs'); hold on;
                xlabel('Normalized $x_1$');
                xlabel('Normalized $x_2$');
            else
                xtx = 1:size(xt1, 2);
                cm = plasma(ceil(1.2*(size(xt1, 1) + size(xt2, 1))));
                hold off;
                for idx = 1:size(xt1, 1)
                    xty = xt1(idx, :);
                    plot(xtx, xty, '-', 'Color', cm(idx, :)); hold on;
                end
                xlabel('Variable number ($\#$)');
                ylabel('Normalized value ($x_\#$)');
            end
            if (size(xt2, 2) == 1)
                xtx = 1:size(xt2, 1);
                plot(xtx, reshape(xt2, 1, numel(xt2)), 'bd'); hold on;
            elseif (size(xt2, 2) == 2)
                plot(xt2(:,1), xt2(:,2), 'bd'); hold on;
            else
                xtx = 1:size(xt2, 2);
                for idx = 1:size(xt2, 1)
                    xty = xt2(idx, :);
                    plot(xtx, xty, '-', 'Color', cm(end-idx+1, :)); hold on;
                end
            end
            if (size(xP, 2) == 1)
                xtx = 1:size(xP, 1);
                plot(xtx, reshape(xP, 1, numel(xP)), 'gh'); hold on;
            elseif (size(xP, 2) == 2)
                plot(xP(:,1), xP(:,2), 'gh'); hold on;
            else
                xtx = 1:size(xP, 2);
                for idx = 1:size(xP, 1)
                    xty = xP(idx, :);
                    plot(xtx, xty, 'g-'); hold on;
                end
            end
            % Objective function space
            if problem.functions.hifi_expensive, return; end
            xt = xt1;
            xt = varScale(xt, problem.bound.xlb, problem.bound.xub, 'unscale');
            [hffF, hffC, hffCEQ] = evaluateHff(problem, xt);
            [idxValid, ~] = separateNaN(xt, hffF, hffC, hffCEQ);
            xt = xt(idxValid, :); hffF = hffF(idxValid, :);
            if (size(hffC, 2) ~= 0), hffC = hffC(idxValid, :);
            else hffC = zeros(size(hffF, 1), 1); end
            hffC = max(hffC, [], 2);
            if (size(hffCEQ, 2) ~= 0), hffCEQ = hffCEQ(idxValid, :);
            else hffCEQ = zeros(size(hffCEQ, 1), 1); end
            hffCEQ = sqrt(max(hffCEQ.^2, [], 2));
            ioptm = (((hffC - problem.control.tolC) <= 0) ...
                & ((hffCEQ - problem.control.tolCEQ) <= 0));
            xts = xt(ioptm, :); xfs = hffF(ioptm, :);
            [xts, xfs, xis] = ndSort(xts, xfs);
            fg_debug2 = debug_figure_open(fg_debug2);
            if (size(xfs, 2) == 1)
                xfx = 1:size(xfs, 1)
                hold off;
                plot(xfx, reshape(xfs, 1, numel(xfs)), 'rs');
                xlabel('Training sample number');
                ylabel('Normalized $f$');
            elseif (size(xfs, 2) == 2)
                hold off;
                plot(xfs(:,1), xfs(:,2), 'rs'); hold on;
                xfss = sortrows(xfs(xis == 1, :), 1);
                plot(xfss(:, 1), xfss(:, 2), 'r-'); hold on;
                xlabel('Normalized $f_1$');
                xlabel('Normalized $f_2$');
            else
                xfx = 1:size(xfs, 2);
                cm = plasma(max(ceil(1.1*size([xt1;xt2], 1)), 1));
                hold off;
                for idx = 1:size(xfs, 1)
                    xfy = xfs(idx, :);
                    plot(xfx, xfy, '-', 'Color', cm(idx, :)); hold on;
                end
                xlabel('Variable number ($\#$)');
                ylabel('Normalized value ($f_\#$)');
            end
            xt = xt2;
            xt = varScale(xt, problem.bound.xlb, problem.bound.xub, 'unscale');
            [hffF, hffC, hffCEQ] = evaluateHff(problem, xt);
            [idxValid, ~] = separateNaN(xt, hffF, hffC, hffCEQ);
            xt = xt(idxValid, :); hffF = hffF(idxValid, :);
            if (size(hffC, 2) ~= 0), hffC = hffC(idxValid, :);
            else hffC = zeros(size(hffF, 1), 1); end
            hffC = max(hffC, [], 2);
            if (size(hffCEQ, 2) ~= 0), hffCEQ = hffCEQ(idxValid, :);
            else hffCEQ = zeros(size(hffCEQ, 1), 1); end
            hffCEQ = sqrt(max(hffCEQ.^2, [], 2));
            ioptm = (((hffC - problem.control.tolC) <= 0) ...
                & ((hffCEQ - problem.control.tolCEQ) <= 0));
            xts = xt(ioptm, :); xfs = hffF(ioptm, :);
            [xts, xfs, xis] = ndSort(xts, xfs);
            fg_debug2 = debug_figure_open(fg_debug2);
            if (size(xfs, 2) == 1)
                xfx = 1:size(xfs, 1)
                plot(xfx, reshape(xfs, 1, numel(xfs)), 'bd');
            elseif (size(xfs, 2) == 2)
                plot(xfs(:,1), xfs(:,2), 'bd'); hold on;
                xfss = sortrows(xfs(xis == 1, :), 1);
                plot(xfss(:, 1), xfss(:, 2), 'b-'); hold on;
            else
                xfx = 1:size(xfs, 2);
                for idx = 1:size(xfs, 1)
                    xfy = xfs(idx, :);
                    plot(xfx, xfy, '-', 'Color', cm(end-idx+1, :)); hold on;
                end
            end
        otherwise
        end
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function fg = debug_figure_open(fg)
    try
        figure(fg);
    catch
        fg = figure('Color',[1 1 1]);
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
