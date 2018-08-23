%% MO-ASMO-II :: ChakongHaimes1983 -- Test problem proposed in Chakong-Haimes 1983
%
% V. Chankong and Y. Y. Haimes,
% "Multiobjective Decision Making Theory and Methodology", Elsevier Science, New York, 1983
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function ChakongHaimes1983_ExpNonlcon
    close all;
    
    % Problem setup
    problem.functions.hifi_obj_exp = @hff_obj;
    problem.functions.hifi_nonlcon_exp = @hff_nonlcon_exp;
    problem.functions.hifi_expensive = false;
    problem.functions.hifi_parallel = true; % if parallel pool > 1, run in parallel
    problem.functions.hifi_vectorized = true; % if no parallel pool, run in vectorized way
    problem.bound.num_x = 2;
    problem.bound.num_f = 2;
    problem.bound.xlb = [-20 -20];
    problem.bound.xub = [20 20];
    problem.bound.flb = [0 -250];
    problem.bound.fub = [250 0];
    problem.sampling.initial.number = 10;
    problem.sampling.update.explore.number = 5;
    problem.sampling.update.exploit.number = 5;
    problem.sampling.validate.number = 5;
    problem.surrogate.method = 'GPR';
    problem.optimization.nsga2.paretofrac = 0.2;
    problem.control.verbose = 2;
    
    % Run MO-ASMO
    problem.control.casefile = mfilename('fullpath');
    resultMOASMO = runMOASMO(problem);
    
    % Run NSGA-II
    problem = resultMOASMO.problem;
    initpop.x = [resultMOASMO.data.c07_poolX_valid{1,1};
                 resultMOASMO.data.c27_valX_valid{1,1}];
    initpop.f = [resultMOASMO.data.c08_poolHffF_valid{1,1};
                 resultMOASMO.data.c29_valHffF_valid{1,1}];
    resultNSGA2 = runDO(problem, 'NSGA-II', initpop);
    
    % Run Epsilon-Constraints
    problem = resultMOASMO.problem;
    xprev = [resultMOASMO.data.c07_poolX_valid{1,1};
             resultMOASMO.data.c27_valX_valid{1,1}];
    fprev = [resultMOASMO.data.c08_poolHffF_valid{1,1};
             resultMOASMO.data.c29_valHffF_valid{1,1}];
    cprev = [resultMOASMO.data.c08_poolHffC_valid{1,1};
             resultMOASMO.data.c29_valHffC_valid{1,1}];
    ceqprev = [resultMOASMO.data.c08_poolHffCEQ_valid{1,1};
               resultMOASMO.data.c29_valHffCEQ_valid{1,1}];
    resultECs = runDO(problem, 'Epsilon-Constraints', xprev, fprev, cprev, ceqprev);
    
    % Save results
    if (exist(problem.control.solpath) == 0)
        mkdir(problem.control.solpath)
    end
    save(fullfile(problem.control.solpath, [problem.control.case, 'results.mat']), ...
        'resultMOASMO', 'resultNSGA2', 'resultECs', '-v7.3');
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function f = hff_obj(x, p)
    f1 = 2 + (x(:,1) - 2).^2 + (x(:,2) - 1).^2;
    f2 = 9*x(:,1) - (x(:,2) - 1).^2;
    f = [f1, f2];
end

function [c, ceq] = hff_nonlcon_exp(x, p)
    c1 = x(:,1) - 3*x(:,2) + 10;
    c2 = x(:,1).^2 + x(:,2).^2 - 225;
    c = [c1, c2];
    ceq = [];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
