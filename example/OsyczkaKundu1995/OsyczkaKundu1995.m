%% MO-ASMO-II :: OsyczkaKundu1995 -- Test problem proposed in Osyczka-Kundu 1995
%
% A. Osyczka and S. Kundu,
% "A new method to solve generalized multicriteria optimization problems using the simple genetic
% algorithm," Structural optimization, 10(2):94-99, 1995.
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function OsyczkaKundu1995
    close all;
    
    % Problem setup
    problem.functions.hifi_obj_exp = @hff_obj;
    problem.functions.hifi_nonlcon_cheap = @hff_nonlcon_cheap;
    problem.functions.hifi_expensive = false;
    problem.functions.hifi_parallel = true; % if parallel pool > 1, run in parallel
    problem.functions.hifi_vectorized = true; % if no parallel pool, run in vectorized way
    problem.bound.num_x = 6;
    problem.bound.num_f = 2;
    problem.bound.xlb = [0 0 1 0 1 0];
    problem.bound.xub = [10 10 5 6 5 10];
    problem.bound.flb = [-300 0];
    problem.bound.fub = [0 100];
    problem.lincon.A = [-1 -1 0 0 0 0;
                        1 1 0 0 0 0;
                        -1 1 0 0 0 0;
                        1 -3 0 0 0 0];
    problem.lincon.b = [-2; 6; 2; 2];
    problem.sampling.initial.number = 12;
    problem.sampling.update.explore.number = 4;
    problem.sampling.update.exploit.number = 2;
    problem.sampling.validation.number = 6;
    problem.stop.residual.ED_avg = 1e-3;
    problem.stop.residual.ED_max = 1e-3;
    problem.stop.residual.HV_data = 'predicted';
    problem.stop.residual.HV_size = 1e-2;
    problem.stop.residual.satisfaction_continuous = 1;
    problem.surrogate.method = 'GPR';
    problem.optimization.nsga2.paretofrac = 0.2;
    problem.optimization.nsga2.functiontolerance = 1e-4;
    problem.optimization.solver = 'NSGA-II';
%     problem.optimization.solver = 'Epsilon-Constraints';
%     problem.optimization.EC.num_per_dim = 30;
%     problem.optimization.EC.obj_num = 2;
%     problem.optimization.fmincon.solver = 'sqp';
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
    f = [-25*(x(:,1)-2).^2 - (x(:,2)-2).^2 - (x(:,3)-1).^2 - (x(:,4)-4).^2 - (x(:,5)-1).^2, ...
        x(:,1).^2 + x(:,2).^2 + x(:,3).^2 + x(:,4).^2 + x(:,5).^2 + x(:,6).^2];
end

function [c, ceq] = hff_nonlcon_cheap(x, p)
    c = [(x(:,3)-3).^2 + x(:,4)-4, ...
        -(x(:,5)-3).^2 - x(:,6)+4];
    ceq = [];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
