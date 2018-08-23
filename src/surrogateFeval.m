%% MO-ASMO-II :: surrogateFeval function
% 1. Evaluate surrogate model for F (and msevalue if DACE is used)
% Usage:
%  [f, msevalue] = surrogateFeval(x, fmodel)
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [f, msevalue] = surrogateFeval(x, fmodel)
    f = [];
    msevalue = [];

    if size(fmodel, 1) == 0
        error('No surrogate model created for objective function');
    end

    method = fmodel.method;
    number = size(x, 1);
    scale_x = fmodel.scale.scale_x;
    scale_f = fmodel.scale.scale_f;

    % Scale input X
    if scale_x
        xlb = fmodel.scale.xlb;
        xub = fmodel.scale.xub;
        x = varScale(x, xlb, xub, 'scale');
    end

    % Evaluate surrogate model for F
    switch lower(method)
    case 'gpr' % Gaussian process regression (Kriging)
        mdl = fmodel.gpr.model;
        num_f = size(mdl, 2);
        f = zeros(number, num_f);
        for idx = 1:num_f
            f(:, idx) = reshape(predict(mdl{idx}, x), number, 1);
        end
    case 'rbf' % Radial-basis function (MO-ASMO-II)
        basisfn = fmodel.rbf.basisfn;
        epsilon = fmodel.rbf.epsilon;
        w = fmodel.rbf.w;
        c = fmodel.rbf.c;
        nc = size(c, 1);
        num_f = size(w, 2);
        f = zeros(number, num_f);
        for idx1 = 1:num_f
            phi = zeros(number, nc);
            for idx2 = 1:nc
                r = sqrt(sum((repmat(c(idx2, :), number, 1) - x).^2, 2));
                switch lower(basisfn)
                case 'linear' % Linear basis fn
                    phi(:,idx2) = r;
                case 'cubic' % Cubic basis fn
                    phi(:,idx2) = r.^3;
                case 'tps' % Thin plate spline basis fn
                    phi(:,idx2) = r.^2.*log(r);
                    phi(r<eps,idx2) = 0;
                case 'gaussian' % Gaussian basis fn
                    phi(:,idx2) = exp(-(epsilon.*r).^2);
                case 'mq' % Multiquadric basis fn
                    phi(:,idx2) = sqrt(1+(epsilon.*r).^2);
                case 'invmq' % Inverse multiquadric basis fn
                    phi(:,idx2) = 1./sqrt(1+(epsilon.*r).^2);
                otherwise
                    error([basisfn,' not supported']);
                end
            end
            f(:, idx1) = phi*w(:, idx1);
        end
    case 'rbn' % Radian-basis function network (MATLAB intrinsic)
        mdl = fmodel.rbn.model;
        f = sim(mdl, x')';
    case 'snn' % Shallow neural network
        mdl = fmodel.snn.model;
        f = sim(mdl, x')';
    case 'dace' % DACE toolkit (external Kriging toolbox)
        mdl = fmodel.dace.model;
        num_f = size(mdl, 2);
        f = zeros(number, num_f);
        msevalue = zeros(number, num_f);
        for idx = 1:num_f
            if size(x, 1) > 1
                [fpredict, msepredict] = predictor(x, mdl{idx});
                f(:, idx) = fpredict(:);
                msevalue(:, idx) = msepredict(:);
            elseif number == 1
                [fpredict, msepredict] = predictor([x; x], mdl{idx});
                f(:, idx) = fpredict(1);
                msevalue(:, idx) = msepredict(1);
            end
        end
    otherwise
        error([method, ' not supported']);
    end

    % Unscale output F
    if scale_f
        flb = fmodel.scale.flb;
        fub = fmodel.scale.fub;
        f = varScale(f, flb, fub, 'unscale');
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
