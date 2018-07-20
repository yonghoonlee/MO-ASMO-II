%% MO-ASMO-II :: trainSurrogate function
% 1. Run high fidelity function
% 2. Obtain objective function and constraint function values
% Usage:
%  surrogate = trainSurrogate(problem, k, xsmp, fsmp)
% Arguments:
%  {problem, x}
% Return:
%  {f, c, ceq}
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function surrogate = trainSurrogate(problem, k, xsmp, fsmp)
    declareGlobalVariables;
    if verbose, disp('Train surrogate models...'); end

    method = problem.surrogate.method;
    number = size(xsmp, 1);
    num_f = problem.bound.num_f;
    xlb = problem.bound.xlb;
    xub = problem.bound.xub;
    flb = problem.bound.flb;
    fub = problem.bound.fub;

    % Adaptive bound adjustment for f (response)
    if (problem.bound.adaptive && (k > 1) && (number > 5))
        flb = min(fsmp, [], 1);
        fub = max(fsmp, [], 1);
    end

    scale_x = problem.surrogate.scale;
    if (scale_x && (min(fub - flb) > 10*eps)), scale_f = true; else, scale_f = false; end

    % Scale x (sample)
    if scale_x
        x = (xsmp - repmat(reshape(xlb, 1, numel(xlb)), number, 1)) ...
            ./ (repmat(reshape(xub, 1, numel(xub)), number, 1) ...
                - repmat(reshape(xlb, 1, numel(xlb)), number, 1));
    else
        x = xsmp;
    end
    % Scale f (response)
    if scale_f
        f = (fsmp - repmat(reshape(flb, 1, numel(flb)), number, 1)) ...
            ./ (repmat(reshape(fub, 1, numel(fub)), number, 1) ...
                - repmat(reshape(flb, 1, numel(flb)), number, 1));
    else
        f = fsmp;
    end

    % Save method and scaling information
    surrogate = [];
    surrogate.method = method;
    surrogate.scale.scale_x = scale_x;
    surrogate.scale.scale_f = scale_f;
    surrogate.scale.xlb = xlb;
    surrogate.scale.xub = xub;
    surrogate.scale.flb = flb;
    surrogate.scale.fub = fub;

    switch lower(method)
    case 'gpr'
    case 'rbf' % Radial-basis function
        basisfn = problem.surrogate.rbf.basisfn;
        epsilon = problem.surrogate.rbf.epsilon;
        c = x;
        w = zeros(number, num_f);
        for k = 1:num_f
            phi = zeros(number, number);
            for j = 1:number
                r = sqrt(sum((repmat(x(j,:), number, 1) - x).^2, 2));
                switch lower(basisfn)
                case 'linear'
                    phi(:, j) = r;
                case 'cubic'
                    phi(:, j) = r.^3;
                case 'tps'
                    phi(:, j) = r.^2.*log(r);
                    phi(r<eps, j) = 0;
                case 'gaussian'
                    phi(:, j) = exp(-(epsilon.*r).^2);
                case 'mq'
                    phi(:, j) = sqrt(1 + (epsilon.*r).^2);
                case 'invmq'
                    phi(:, j) = 1./sqrt(1 + (epsilon.*r).^2);
                otherwise
                    error([basisfn, ' not supported']);
                end
            end
            w(:, k) = pinv(phi)*f(:, k);
        end
        surrogate.rbf.basisfn = basisfn;
        surrogate.rbf.epsilon = epsilon;
        surrogate.rbf.w = w;
        surrogate.rbf.c = c;
    case 'rbn'
    case 'snn'
    case 'dace'
    otherwise
    end

    




end