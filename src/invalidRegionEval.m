%% MO-ASMO-II :: invalidRegionEval function
% 1. Evaluate invalid region model for given input
% Usage:
%  c_invalid = invalidRegionEval(x, irmodel)
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function c_invalid = invalidRegionEval(x, irmodel)
    % Initially return variable has a fixed size output (to prevent error in optimization)
    c_invalid = zeros(size(x, 1), 1);

    % Check if irmodel exists
    if (size(irmodel, 1) == 0), return; end

    % Scale input X
    scale_x = irmodel.scale.scale_x;
    xlb = irmodel.scale.xlb;
    xub = irmodel.scale.xub;
    if scale_x, x = varScale(x, xlb, xub, 'scale'); end

    p = irmodel.parameter;
    number = size(p.x, 1);
    numberx = size(x, 1);
    dopt = irmodel.dopt;

    for idx = 1:numberx
        z = x(idx, :);
        Kxz = sum(dopt.*exp(-p.q*sum((p.x - repmat(z, number, 1)).^2, 2)));
        Kxx = sum(dopt.*exp(-p.q*sum((p.x - repmat(p.x(1, :), number, 1)).^2, 2)));
        c_invalid(idx, 1) = 2*(Kxz - Kxx) + p.epsilon;
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
