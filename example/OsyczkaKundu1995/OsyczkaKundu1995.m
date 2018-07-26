%% MO-ASMO-II :: OsyczkaKundu1995 -- Test problem proposed in Osyczka-Kundu 1995
%
% A. Osyczka and S. Kundu,
% "A new method to solve generalized multicriteria optimization problems using the simple genetic
% algorithm," Structural optimization, 10(2):94-99, 1995.
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function OsyczkaKundu1995
    % Problem setup
    problem.functions.hifi_obj_exp = @hff_obj;
    problem.functions.hifi_nonlcon_cheap = @hff_nonlcon_cheap;
    problem.functions.hifi_expensive = false;
    problem.functions.hifi_parallel = true; % if parallel pool > 1, run in parallel
    problem.functions.hifi_vectorized = true; % if no parallel pool, run in vectorized way
    problem.bound.num_x = 6 ;
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
    problem.sampling.initial.number = 10;
    problem.surrogate.method = 'GPR';
    % Run MO-ASMO
    problem.control.casefile = mfilename('fullpath');
    result1 = runMOASMO(problem);
    % Run NSGA-II
    problem = result1.problem;
    initpop.x = [result1.data.c07_poolX_valid{1,1}; result1.data.c27_valX_valid{1,1}];
    initpop.f = [result1.data.c08_poolHffF_valid{1,1}; result1.data.c29_valHffF_valid{1,1}];
    result2 = runDO(problem, 'NSGA-II', initpop);
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
