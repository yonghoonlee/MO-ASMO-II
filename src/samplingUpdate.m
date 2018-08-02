%% MO-ASMO-II :: samplingUpdate function
% 1. Generate update samples for updating surrogate model
% Usage:
%  xt = samplingUpdate(varargin)
% Arguments:
%  {problem, k, poolXvalid, irmodel, xP}
%
% Force-directed layout algorithm used in exploitation sampling of predicted Pareto set. 
% Algorithm came from a combination and modification of following references:
% [1] Peter Eades, "A heuristic for graph drawing," Congressus Numerantium, 42:149-160, 1984.
% [2] Thomas M. J. Fruchterman and Edward M. Reingold, "Graph drawing by force-directed placement,"
%     Software: Practice and Experience, 21(11):1129-1164, 1991.
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function xt = samplingUpdate(problem, k, poolXvalid, poolXinvalid, xP)
    declareGlobalVariables;

    exploit_method = problem.sampling.update.exploit.method;
    exploit_number = problem.sampling.update.exploit.number;
    explore_method = problem.sampling.update.explore.method;
    explore_number = problem.sampling.update.explore.number;
    dimension = problem.bound.num_x;
    xlb = problem.bound.xlb;
    xub = problem.bound.xub;

    % Exploration
    if (explore_number > 0)
        ex = -1; % Run exploration sampling
    else
        ex = 2; % Bypass exploration sampling
        xt1 = [];
    end
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
        irmodel = trainInvalidRegion(problem, poolXinvalid);
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
                [value, flg] = fmincon_explore(problem, xt1(idx, :), poolXvalid, irmodel);
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
    if (exploit_number > 0)
        if verbose, disp('Exploitation sampling...'); end
        switch lower(exploit_method)
        case 'fdl'
            % Force-driected layout method
            number_xP = size(xP, 1);
            num_xP = size(xP, 2);
            % Scale
            xPlb = min(xP,[],1);
            xPub = max(xP,[],1); % set bounding box for Pareto set
            xPlb = xPlb - 0.1*(xPub - xPlb);
            xPub = xPub + 0.1*(xPub - xPlb); % expand 10% below and above bounds
            if sum(xPub - xPlb) == 0 % If ub and lb are the same
                xPub = xPub + 0.1*(reshape(xub, 1, numel(xub)) - reshape(xlb, 1, numel(xlb)));
                xPlb = xPlb - 0.1*(reshape(xub, 1, numel(xub)) - reshape(xlb, 1, numel(xlb)));
            end
            xPlb = max([xPlb; reshape(xlb, 1, numel(xlb))]);
            xPub = min([xPub; reshape(xub, 1, numel(xub))]); % comply original bound
            xPS = varScale(xP, xPlb, xPub, 'scale');
            % Find cluster centers
            if (exploit_number >= number_xP), exploit_number = number_xP - 1; end
            if (exploit_number <= 1)
                xB = xPS;
            else
                [~, xB, ~] = kmeans(xPS, exploit_number);
            end
            number_B = size(xB, 1);
            % Generate nearby points to the xB
            pm = rand(size(xB));
            pm(pm<0.5) = -1; pm(pm>=0.5) = 1;
            xt2 = xB + 0.05*pm;
            % Main loop of force-directed layout algorithm
            k = 1/number_B;
            Rmax = 1e-9;
            imax = 10;
            for iloop = 1:imax
                tcool = imax/iloop;
                % Run FDL
                Vdisp = zeros(size(xt2)); % Displacement of vertices
                % Attractive and repulsive forces between Base point (xB) and Moving point (xM)
                Xdiff = xt2 - xB;
                Rdiff = sqrt(sum(Xdiff.^2, 2));
                Vdisp = Vdisp - (Xdiff./Rdiff).*log10(100*Rdiff);
                % Repulsive forces between Pareto set (xPS) and Moving point (xM)
                for idx = 1:number_xP
                    Xdiff = xt2 - repmat(xPS(idx, :), number_B, 1);
                    Rdiff = sqrt(sum(Xdiff.^2, 2));
                    Vdisp = Vdisp + (Xdiff./Rdiff).*(k^2./(Rdiff*10).^2);
                end
                % Limit maximum displacement
                R = sqrt(sum(Vdisp.^2, 2));
                Rmax = max(Rmax, max(R));
                xt2 = xt2 + Vdisp./R.*min([R, repmat(tcool, number_B, 1)], [], 2)/(Rmax*10);
                for idx = 1:size(xt2, 2)
                    VC = xt2(:, idx);
                    VC(VC<0) = 0;
                    VC(VC>1) = 1;
                    xt2(:,idx) = VC;
                end
                if (verbose == 2), debugAnalysis(problem, 'force-directed', xt2, xB, xPS); end
            end
            % Unscale
            xt2 = varScale(xt2, xPlb, xPub, 'unscale');
        case 'cdd'
        end
    else
        xt2 = [];
    end

    % Combine
    xt = [xt1; xt2];
    xt = unique(xt, 'rows');

    if (verbose == 2), debugAnalysis(problem, 'UpdateSampling', xt1, xt2, xP); end
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
        f2 = sum(-log(max(sum((x - xprev).^2, 2), 1e-12))) / size(xprev, 1);
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
