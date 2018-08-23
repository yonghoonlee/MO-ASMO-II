%% MO-ASMO-II :: Lee2017a -- Test problem proposed in Lee et al., IDETC 2017
%
% Yong Hoon Lee, R. E. Corman, Randy H. Ewoldt, and James T. Allison,
% 
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function Lee2017a_NonlconExp
    close all;
    
    % Problem setup
    problem.functions.hifi_obj_exp = @hff_obj;
    problem.functions.hifi_nonlcon_exp = @hff_nonlcon_exp;
    problem.functions.hifi_expensive = false;
    problem.functions.hifi_parallel = true; % if parallel pool > 1, run in parallel
    problem.functions.hifi_vectorized = true; % if no parallel pool, run in vectorized way
    problem.bound.num_x = 2;
    problem.bound.num_f = 2;
    problem.bound.xlb = [-5 -5];
    problem.bound.xub = [5 5];
    problem.bound.flb = [-1 -1];
    problem.bound.fub = [7 0];
    problem.sampling.initial.number = 30;
    problem.sampling.update.explore.number = 10;
    problem.sampling.update.explore.w1_distance = 1e-3;
    problem.sampling.update.explore.w2_disperse = 1e-9;
    problem.sampling.update.exploit.number = 10;
    problem.sampling.validate.number = 10;
    problem.surrogate.method = 'GPR';
    problem.stop.maxiter = 100;
    problem.stop.residual.satisfaction_continuous = 5;
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
    
    % Save results
    if (exist(problem.control.solpath) == 0)
        mkdir(problem.control.solpath)
    end
    save(fullfile(problem.control.solpath, [problem.control.case, 'results.mat']), ...
        'resultMOASMO', 'resultNSGA2', '-v7.3');
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function f = hff_obj(x, p)
    f1 = (3*sin(2.5*x(:,1)) - 2*x(:,1)).*cos(x(:,2)).*exp(-(1e-3)*x(:,2).^2);
    f2 = (3/8)*(0.3*abs(x(:,1))).^(19/25).*(sin(5*x(:,2)).*x(:,2));
    f = [f1, f2];
end

function [c, ceq] = hff_nonlcon_exp(x, p)
    c = (100*(x(:,2) - x(:,1).^2).^2+(x(:,1) - 1).^2) - 1;
    ceq = [];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
