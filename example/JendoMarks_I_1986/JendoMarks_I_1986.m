%% MO-ASMO-II :: JendoMarks_I_1986 -- Test problem proposed in section 4.1 of Jendo-Marks 1986
% S. Jendo and W. Marks,
% "Multiobjective structural optimization,"
% In Prékopa A., Szelezsáan J., Strazicky B. (eds)
% System Modelling and Optimization. Lecture Notes in Control and Information Sciences,
% vol 84. Springer, Berlin, Heidelberg
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function JendoMarks_I_1986
    close all;
    
    % Problem setup
    problem.functions.hifi_combined_exp = @hff_combined;
    problem.functions.hifi_expensive = false;
    problem.functions.hifi_parallel = true; % if parallel pool > 1, run in parallel
    problem.functions.hifi_vectorized = false; % if no parallel pool, run in vectorized way
    problem.parameter.sigma1.lb = -500;
    problem.parameter.sigma1.ub = 10000;
    problem.parameter.sigma2.lb = -500;
    problem.parameter.sigma2.ub = 10000;
    problem.parameter.M1 = 320;
    problem.parameter.M2 = 100;
    problem.parameter.d = 0.07;
    problem.parameter.k = 0.000054;
    problem.bound.num_x = 6;
    problem.bound.num_f = 2;
    problem.bound.xlb = [0.30 0.10 0.10 0.08 100 0.08];
    problem.bound.xub = [0.70 0.40 0.20 0.24 1000 0.28];
    problem.bound.flb = [0.10 600];
    problem.bound.fub = [0.15 800];
    problem.sampling.initial.number = 15;
    problem.sampling.update.explore.number = 10;
    problem.sampling.update.exploit.number = 5;
    problem.sampling.validate.number = 5;
    problem.surrogate.method = 'GPR';
    problem.optimization.nsga2.paretofrac = 0.2;
    problem.optimization.nsga2.maxgen = 200;
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

function [f,c,ceq] = hff_combined(x, p)
    % Design variables
    x1 = x(1); % beam depth
    x2 = x(2); % flange width
    x3 = x(3); % flange thickness
    x4 = x(4); % web thickness
    x5 = x(5); % prestressing force
    x6 = x(6); % eccentricity of the prestressing force
    Jf = (1/12)*x2*x3^3; % flange moment of inertia
    Af = x2*x3;
    df = x1/2 - x3/2; % flange distance
    Jw = (1/12)*x4*(x1 - 2*x3)^3; % web moment of inertia
    Aw = (x1 - 2*x3)*x4;
    dw = 0;
    J = 2*(Jf + Af*df^2) + (Jw + Aw*dw^2);
    % Minimum volume of the concrete
    A = 2*x2*x3 + x4*(x1 - 2*x3);
    % Prestressing steel
    S = x5;
    C1 = x5/A - x5*x6*x1/(2*J) + p.M1*x1/(2*J) - p.sigma1.ub;
    C2 = -(x5/A + x5*x6*x1/(2*J) - p.M1*x1/(2*J) - p.sigma1.lb);
    C3 = -(x5/A - x5*x6*x1/(2*J) + p.M2*x1/(2*J) - p.sigma2.lb);
    C4 = x5/A + x5*x6*x1/(2*J) - p.M2*x1/(2*J) - p.sigma2.ub;
    C5 = x6 - x1/2 + p.d;
    f = [A, S];
    c = [C1, C2, C3, C4, C5];
    ceq = [];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
