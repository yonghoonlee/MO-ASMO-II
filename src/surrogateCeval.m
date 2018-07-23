%% MO-ASMO-II :: surrogateCeval function
% 1. Evaluate surrogate model for C and CEQ (and msevalue if DACE is used)
% Usage:
%  [c, ceq, c_msevalue, ceq_msevalue] = surrogateCeval(x, cmodel, ceqmodel)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [c, ceq, c_msevalue, ceq_msevalue] = surrogateCeval(x, cmodel, ceqmodel, irmodel)
    c = [];
    ceq = []
    c_invalid = [];
    c_msevalue = [];
    ceq_msevalue = [];

    if ((size(cmodel, 1) == 0) && (size(ceqmodel, 1) == 0)), return; end
    
    % Evaluate C
    xinput = x; % backup input X for CEQ evaluation later
    if (size(cmodel, 1) ~= 0)
        method = cmodel.method;
        number = size(x, 1);
        scale_x = cmodel.scale.scale_x;
        scale_f = cmodel.scale.scale_f;

        % Scale input X
        if scale_x
            xlb = cmodel.scale.xlb;
            xub = cmodel.scale.xub;
            x = varScale(x, xlb, xub, 'scale');
        end

        % Evaluate surrogate model for C
        switch lower(method)
        case 'gpr' % Gaussian process regression (Kriging)
            mdl = cmodel.gpr.model;
            num_c = size(mdl, 2);
            c = zeros(number, num_c);
            for idx = 1:num_c
                c(:, idx) = reshape(predict(mdl{idx}, x), number, 1);
            end
        case 'rbf' % Radial-basis function (MO-ASMO-II)
            basisfn = cmodel.rbf.basisfn;
            epsilon = cmodel.rbf.epsilon;
            w = cmodel.rbf.w;
            ctr = cmodel.rbf.c;
            nctr = size(ctr, 1);
            num_c = size(w, 2);
            c = zeros(number, num_c);
            for idx1 = 1:num_c
                phi = zeros(number, nctr);
                for idx2 = 1:nctr
                    r = sqrt(sum((repmat(ctr(idx2, :), number, 1) - x).^2, 2));
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
                c(:, idx1) = phi*w(:, idx1);
            end
        case 'rbn' % Radian-basis function network (MATLAB intrinsic)
            mdl = cmodel.rbn.model;
            c = sim(mdl, x')';
        case 'snn' % Shallow neural network
            mdl = cmodel.snn.model;
            c = sim(mdl, x')';
        case 'dace' % DACE toolkit (external Kriging toolbox)
            mdl = cmodel.dace.model;
            num_c = size(mdl, 2);
            c = zeros(number, num_c);
            msevalue = zeros(number, num_c);
            for idx = 1:num_c
                if size(x, 1) > 1
                    [c_predict, c_msepredict] = predictor(x, mdl{idx});
                    c(:, idx) = c_predict(:);
                    c_msevalue(:, idx) = c_msepredict(:);
                elseif number == 1
                    [c_predict, c_msepredict] = predictor([x; x], mdl{idx});
                    c(:, idx) = c_predict(1);
                    c_msevalue(:, idx) = c_msepredict(1);
                end
            end
        otherwise
            error([method, ' not supported']);
        end

        % Descale output C
        if scale_f
            flb = cmodel.scale.flb;
            fub = cmodel.scale.fub;
            c = varScale(c, flb, fub, 'descale');
        end
    end

    % Evaluate CEQ
    x = xinput; % restore input X and discarding previously scaled one
    if (size(ceqmodel, 1) ~= 0)
        method = ceqmodel.method;
        number = size(x, 1);
        scale_x = ceqmodel.scale.scale_x;
        scale_f = ceqmodel.scale.scale_f;

        % Scale input X
        if scale_x
            xlb = ceqmodel.scale.xlb;
            xub = ceqmodel.scale.xub;
            x = varScale(x, xlb, xub, 'scale');
        end

        % Evaluate surrogate model for CEQ
        switch lower(method)
        case 'gpr' % Gaussian process regression (Kriging)
            mdl = ceqmodel.gpr.model;
            num_ceq = size(mdl, 2);
            ceq = zeros(number, num_ceq);
            for idx = 1:num_ceq
                ceq(:, idx) = reshape(predict(mdl{idx}, x), number, 1);
            end
        case 'rbf' % Radial-basis function (MO-ASMO-II)
            basisfn = ceqmodel.rbf.basisfn;
            epsilon = ceqmodel.rbf.epsilon;
            w = ceqmodel.rbf.w;
            ctr = ceqmodel.rbf.c;
            nctr = size(ctr, 1);
            num_ceq = size(w, 2);
            ceq = zeros(number, num_ceq);
            for idx1 = 1:num_ceq
                phi = zeros(number, nctr);
                for idx2 = 1:nctr
                    r = sqrt(sum((repmat(ctr(idx2, :), number, 1) - x).^2, 2));
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
                ceq(:, idx1) = phi*w(:, idx1);
            end
        case 'rbn' % Radian-basis function network (MATLAB intrinsic)
            mdl = ceqmodel.rbn.model;
            ceq = sim(mdl, x')';
        case 'snn' % Shallow neural network
            mdl = ceqmodel.snn.model;
            ceq = sim(mdl, x')';
        case 'dace' % DACE toolkit (external Kriging toolbox)
            mdl = ceqmodel.dace.model;
            num_ceq = size(mdl, 2);
            ceq = zeros(number, num_ceq);
            msevalue = zeros(number, num_ceq);
            for idx = 1:num_ceq
                if size(x, 1) > 1
                    [ceq_predict, ceq_msepredict] = predictor(x, mdl{idx});
                    ceq(:, idx) = ceq_predict(:);
                    ceq_msevalue(:, idx) = ceq_msepredict(:);
                elseif number == 1
                    [ceq_predict, ceq_msepredict] = predictor([x; x], mdl{idx});
                    ceq(:, idx) = ceq_predict(1);
                    ceq_msevalue(:, idx) = ceq_msepredict(1);
                end
            end
        otherwise
            error([method, ' not supported']);
        end

        % Descale output CEQ
        if scale_f
            flb = ceqmodel.scale.flb;
            fub = ceqmodel.scale.fub;
            ceq = varScale(ceq, flb, fub, 'descale');
        end
    end

    % Evaluate invalid input zone constraints
    x = xinput;
    c_invalid = invalidRegionEval(x, irmodel);

    % Concatenate C and C_INVALID
    c = [c, c_invalid];
    if (size(c_msevalue, 1) ~= 0)
        c_msevalue = [c_msevalue, zeros(size(c_invalid))];
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
