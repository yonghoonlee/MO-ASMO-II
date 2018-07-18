%% MO-ASMO-II :: parallelPoolSize function
% 1. Returns parallel pool size if exists.
% 2. Otherwise, returns 0.
% Usage:
%  npool = parallelPoolSize()
% Return:
%  {npool}
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function npool = parallelPoolSize()
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    if isempty(poolobj)
        npool = 0;
    else
        npool = poolobj.NumWorkers;
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
