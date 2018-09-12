%% MO-ASMO-II :: enforceIndexLincon function
% 1. Among x data of iHff = 1, make iHff = 0 for those that does not comply linear constraints
% Usage:
%  iHff = enforceIndexLincon(problem, iHff, x)
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function iHff = enforceIndexLincon(problem, iHff, x)
    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    
    iHff = logical(iHff);
    if (size(b, 1) == 0) && (size(beq, 1) == 0)
        return;
    end
    
    for idx = 1:size(x, 1)
        xsel = reshape(x(idx, :), size(x, 2), 1);
        if (size(b, 1) ~= 0)
            C = A*xsel - b - problem.control.tolC;
            if (max(C) <= 0) && (iHff(idx) == true)
                iHff(idx) = true;
            else
                iHff(idx) = false;
            end
        end
        if (size(beq, 1) ~= 0)
            CEQ = sqrt((Aeq*xsel - beq).^2) - problem.control.tolCEQ;
            if (max(CEQ) <= 0) && (iHff(idx) == true)
                iHff(idx) = true;
            else
                iHff(idx) = false;
            end
        end
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
