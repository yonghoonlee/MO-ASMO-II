%% MO-ASMO-II :: runMOASMO function
% 1. Read problem structure
% 2. Run MO-ASMO-II algorithm
% 3. Return result structure
% Usage:
%  result = runMOASMO(varargin)
% Arguments:
%  * Recommended:
%   - {problem}
%   - {result} -- for restart problem. problem structure should be a field in the result structure
%  * Obsolete (but supported):
%   - {(@settings), @hifi_combined, casefile}
%   - {(@settings), @hifi_obj, (@hifi_nonlcon), casefile}
%   - {(@settings), @hifi_combined, casefile, outeriter, innercase}
%   - {(@settings), @hifi_obj, (@hifi_nonlcon), casefile, outeriter, innercase}
% Note:
%  Functions @hifi_obj, @hifi_nonlcon, @hifi_combined are assumed expensive.
%  To include cheap constraints that can be evaluated prior to the objective function evaluation,
%  please use problem structure as the only input argument.
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function result = runMOASMO(varargin)
    declareGlobalVariables;
    
    % Parse input arguments (varargin). Run "help runMOASMO" command for more information
    [problem, result, restart] = parseInputArguments(varargin{:});

    % If result is provided as an input, restart problem
    if restart
        data = result.data;
        scriptExtractDataStructure;
        k = size(data.c37_resHVRhffHistory{1,1}, 1);
    else
        problem = problemSetup(problem);
        result.problem = problem;
        k = 0;
    end
    
    if verbose, disp('Multi-Objective Adaptive Surrogate Model-Based Optimization (MO-ASMO)'); end

    % Main loop for MO-ASMO algorithm
    while (true)
        k = k + 1;
        if verbose, disp(['[[Iteration:',num2str(k),']]']); end
        time_elapsed = [];
        
        % Initial / update sampling
        prepRandomSeed();
        if k == 1
            [c01_smpX, t_elapsed] = samplingInitial(problem);
        else
            [c01_smpX, t_elapsed] = samplingUpdate( ...
                problem, k, c07_poolX_valid, c09_poolX_invalid, c14_parX);
        end
        time_elapsed.c01_sampling = t_elapsed;

        % High fidelity function evaluation
        [c02_hffF, c02_hffC, c02_hffCEQ, t_elapsed] = evaluateHff(problem, c01_smpX);
        [idxValid, idxInvalid] = separateNaN(c01_smpX, c02_hffF, c02_hffC, c02_hffCEQ);
        scriptSeparateSmpHff;
        time_elapsed.c02_evaluateHff = t_elapsed;

        % Merge pervious and current results to generate training point pool
        scriptMergePrevAndCurrentPool;
        
        % Surrogate model construction
        [c11_surrogateF, t_elapsed_F] = trainSurrogate( ...
            problem, c07_poolX_valid, c08_poolHffF_valid, true);
        [c11_surrogateC, t_elapsed_C] = trainSurrogate( ...
            problem, c07_poolX_valid, c08_poolHffC_valid, false);
        [c11_surrogateCEQ, t_elapsed_CEQ] = trainSurrogate( ...
            problem, c07_poolX_valid, c08_poolHffCEQ_valid, false);
        time_elapsed.c07_trainSurrogateF = t_elapsed_F;
        time_elapsed.c07_trainSurrogateC = t_elapsed_C;
        time_elapsed.c07_trainSurrogateCEQ = t_elapsed_CEQ;

        % Invalid region model construction
        if size(c09_poolX_invalid,1) == 0
            c12_invalidRegion = [];
            t_elapsed = 0;
        else
            [c12_invalidRegion, t_elapsed] = trainInvalidRegion(problem, c09_poolX_invalid);
        end
        time_elapsed.c12_trainIRmodel = t_elapsed;

        % Generate starting points for optimization using previous results
        if (k == 1), c14_parX = []; c15_parSurF = []; end
        c13_startingPoints = compileStartingPoints(problem, ...
            c07_poolX_valid, c08_poolHffF_valid, c08_poolHffC_valid, c08_poolHffCEQ_valid, ...
            c14_parX, c15_parSurF);

        % Surrogate model-based optimization
        [c14_parX, c15_parSurF, c15_parSurC, c15_parSurCEQ, c16_parSurOut, t_elapsed] ...
            = surrogateOptim(problem, c11_surrogateF, c11_surrogateC, c11_surrogateCEQ, ...
                c12_invalidRegion, c13_startingPoints);
        time_elapsed.c14_SBO = t_elapsed;

        % (Optional) Compute exact error if high fidelity function evaluation is not expensive
        if ~(problem.functions.hifi_expensive)
            [c17_parHffF, c17_parHffC, c17_parHffCEQ] = evaluateHff(problem, c14_parX);
        else
            c17_parHffF = c15_parSurF;
            c17_parHffC = c15_parSurC;
            c17_parHffCEQ = c15_parSurCEQ;
        end
        [idxValid, idxInvalid] ...
            = separateNaN(c14_parX, c17_parHffF, c17_parHffC, c17_parHffCEQ);
        scriptSeparateParSurHff;

        % Validation sampling, high fidelity function evaluations
        [c24_valX, c25_valSurF, c25_valSurC, c25_valSurCEQ, t_elapsed_smp] ...
            = samplingValidation(problem, c14_parX, c15_parSurF, c15_parSurC, c15_parSurCEQ);
        [c26_valHffF, c26_valHffC, c26_valHffCEQ, t_elapsed_hff] = evaluateHff(problem, c24_valX);
        [idxValid, idxInvalid] ...
            = separateNaN(c24_valX, c26_valHffF, c26_valHffC, c26_valHffCEQ);
        scriptSeparateValSurHff;
        time_elapsed.c24_valSmp = t_elapsed_smp;
        time_elapsed.c26_valHff = t_elapsed_hff;

        % Compute convergence and error metrices
        % Eulerian distance (ED) error for F
        if (k == 1), c33_prevEDarrF = []; else, c33_prevEDarrF{k - 1, 1} = c33_EDvecF; end
        [c33_EDvecF, c33_EDavgF] = errEuclideanDistance(problem, ...
            c29_valHffF_valid, c28_valSurF_valid, c33_prevEDarrF, c11_surrogateF);
        % Hypervolume (HV), hypervolume ratio (HVR), Normalized residuals
        scriptHVResidualParetoPred; % of predicted Pareto set
        % Hypervolume (HV) and hypervolume ratio (HVR), Normalized residuals
        scriptHVResidualParetoHff; % of Pareto set of high-fidelity results

        % Compile data to data structure
        scriptCompileDataStructure;
        result.data = []; result.data = data;
        
        % Compile solution to result structure, count high-fidelity function evaluations
        scriptCompileSolutionToResult;
        
        % Compile other analysis metrics to result structure
        time_elapsed = struct2table(time_elapsed);
        time_elapsed.iter_total = sum(time_elapsed{1,:}, 2);
        result.time.elapsed_iter = [result.time.elapsed_iter; time_elapsed];
        result.time.elapsed_total = time_elapsed;
        result.time.elapsed_total{1,:} = sum(result.time.elapsed_iter{:,:}, 1);
        result.problem.control.randomseed = randomseed;
        
        % Save intermediate result
        scriptSaveIntermediateResult; % save to .mat file

        % Stopping criteria evaluation
        if (exist('stopcount') ~= 1), stopcount = 0; end
        [stoploop, stopcount] = evaluateStopCriteria( ...
            problem, k, stopcount, c33_EDvecF, c33_prevEDarrF, ...
            c35_resHVpred, c35_resHVRpred, c37_resHVhff, c37_resHVRhff, ...
            c35_resHVpredHistory, c35_resHVRpredHistory, c37_resHVhffHistory, c37_resHVRhffHistory);
        if stoploop, break; end
    end

end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
