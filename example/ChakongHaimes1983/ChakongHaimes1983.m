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

function ChakongHaimes1983
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
    problem.sampling.initial.number = 4;
    problem.sampling.update.explore.number = 2;
    problem.sampling.update.exploit.number = 2;
    problem.sampling.validation.number = 4;
    problem.stop.residual.ED_avg = 1e-3;
    problem.stop.residual.ED_max = 1e-3;
    problem.stop.residual.HV_data = 'predicted';
    problem.stop.residual.HV_size = 1e-2;
    problem.stop.residual.satisfaction_continuous = 1;
    problem.surrogate.method = 'GPR';
    problem.optimization.nsga2.paretofrac = 0.2;
    problem.optimization.solver = 'Epsilon-Constraints';
    problem.optimization.EC.num_per_dim = 24;
    problem.optimization.fmincon.solver = 'sqp';
    problem.control.verbose = 2;
    
    % Run MO-ASMO
    problem.control.casefile = mfilename('fullpath');
    resultMOASMO = runMOASMO(problem);
    
    % Save results
    if (exist(problem.control.solpath) == 0)
        mkdir(problem.control.solpath)
    end
    save(fullfile(problem.control.solpath, [problem.control.case, 'results.mat']), ...
        'resultMOASMO', '-v7.3');
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function f = hff_obj(x, p)
    f1 = 2 + (x(:,1) - 2).^2 + (x(:,2) - 1).^2;
    f2 = 9*x(:,1) - (x(:,2) - 1).^2;
    f = [f1, f2];
end

function [c, ceq] = hff_nonlcon_cheap(x, p)
    c = x(:,1).^2 + x(:,2).^2 - 225;
    ceq = [];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
