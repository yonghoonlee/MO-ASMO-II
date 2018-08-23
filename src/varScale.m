%% MO-ASMO-II :: varScale function
% 1. Scale and unscale variables regarding lower (xlb) and upper (xub) bounds.
% Usage:
%  xout = varScale(xin, xlb, xub, method)
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function xout = varScale(xin, xlb, xub, method)
    repeat = size(xin, 1);
    if (repeat == 0), xout = xin; return; end
    xlb = repmat(reshape(xlb, 1, numel(xlb)), repeat, 1);
    xub = repmat(reshape(xub, 1, numel(xub)), repeat, 1);

    switch lower(method)
    case 'scale'
        xout = (xin - xlb)./(xub - xlb);
    case 'unscale'
        xout = xlb + (xub - xlb).*xin;
    otherwise
        error('Function samplingScale requires method of either scale or unscale');
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
