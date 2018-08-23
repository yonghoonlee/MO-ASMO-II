%% MO-ASMO-II :: scriptSaveIntermediateResult script
% 1. Save intermediate data stored in result structure to .mat file
% Usage:
%  scriptSaveIntermediateResult
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

if (problem.nested.outeriter == 0) && (problem.nested.innercase == 0)
    if (exist(problem.control.solpath) == 0)
        mkdir(problem.control.solpath)
    end
    save(fullfile(problem.control.solpath, ...
        [problem.control.case, '_iter', num2str(k,'%04d'), '.mat']), ...
        'result', '-v7.3');
else
    if (exist(problem.control.solpath) == 0)
        mkdir(problem.control.solpath)
    end
    save(fullfile(problem.control.solpath, ...
        [problem.control.case, ...
        '_outeriter', num2str(problem.nested.outeriter,'%04d'), ...
        '_innercase', num2str(problem.nested.innercase,'%04d'), ...
        '_iter', num2str(k,'%04d'), '.mat']), ...
        'result', '-v7.3');
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
