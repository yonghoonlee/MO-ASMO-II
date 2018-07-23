%% MO-ASMO-II :: runMOASMO function
% 1. Read problem structure
% 2. Run MO-ASMO-II algorithm
% 3. Return result structure
% Usage:
%  result = runMOASMO(varargin)
% Arguments:
%  {problem}
%  {(@settings), @hifi_combined, casefile}
%  {(@settings), @hifi_obj, (@hifi_nonlcon), casefile}
%  {(@settings), @hifi_combined, casefile, outeriter, innercase}
%  {(@settings), @hifi_obj, (@hifi_nonlcon), casefile, outeriter, innercase}
% Note:
%  Functions @hifi_obj, @hifi_nonlcon, @hifi_combined are assumed expensive.
%  To include cheap constraints that can be evaluated prior to the objective function evaluation,
%  please use problem structure as the only input argument.
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function result = runMOASMO(varargin)
    declareGlobalVariables;
    
    [problem, result, restart] = parseInputArguments(varargin{:});

    if restart
        data = result.data;
        scriptExtractDataStructure;
    else
        problem = problemSetup(problem);
        result.problem = problem;
        k = 0;
    end

    while (true)
        k = k + 1;
        if verbose, disp(['[[Iteration:',num2str(k),']]']); end
        
        % Initial / update sampling
        prepRandomSeed();
        if k == 1
            c01_smpX = samplingInitial(problem);
        else
            c01_smpX = samplingUpdate(problem, k, c07_poolX_valid, c09_poolX_invalid, c14_parX);
        end

        % High fidelity function evaluation
        [c02_hffF, c02_hffC, c02_hffCEQ] = evaluateHff(problem, c01_smpX);
        [idxValid, idxInvalid] = separateNaN(c01_smpX, c02_hffF, c02_hffC, c02_hffCEQ);
        c03_smpX_valid = c01_smpX(idxValid,:);                  % x         valid
        c04_hffF_valid = c02_hffF(idxValid,:);                  % f         valid
        c04_hffC_valid = c02_hffC(idxValid,:);                  % c         valid
        c04_hffCEQ_valid = c02_hffCEQ(idxValid,:);              % ceq       valid
        c05_smpX_invalid = c01_smpX(idxInvalid,:);              % x         invalid
        c06_hffF_invalid = c02_hffF(idxInvalid,:);              % f         invalid
        c06_hffC_invalid = c02_hffC(idxInvalid,:);              % c         invalid
        c06_hffCEQ_invalid = c02_hffCEQ(idxInvalid,:);          % ceq       invalid

        % Merge pervious and current results to generate training point pool
        if k == 1
            % valid data
            c07_poolX_valid = c03_smpX_valid;
            c08_poolHffF_valid = c04_hffF_valid;
            c08_poolHffC_valid = c04_hffC_valid;
            c08_poolHffCEQ_valid = c04_hffCEQ_valid;
            % invalid data
            c09_poolX_invalid = c05_smpX_invalid;
            c10_poolHffF_invalid = c06_hffF_invalid;
            c10_poolHffC_invalid = c06_hffC_invalid;
            c10_poolHffCEQ_invalid = c06_hffCEQ_invalid;
        else
            % valid data
            c07_poolX_valid = [c07_poolX_valid; c27_valX_valid; c03_smpX_valid];
            c08_poolHffF_valid = [c08_poolHffF_valid; c29_valHffF_valid; c04_hffF_valid];
            c08_poolHffC_valid = [c08_poolHffC_valid; c29_valHffC_valid; c04_hffC_valid];
            c08_poolHffCEQ_valid = [c08_poolHffCEQ_valid; c29_valHffCEQ_valid; c04_hffCEQ_valid];
            % invalid data
            c09_poolX_invalid ...
                = [c09_poolX_invalid; c30_valX_invalid; c05_smpX_invalid];
            c10_poolHffF_invalid ...
                = [c10_poolHffF_invalid; c32_valHffF_invalid; c06_hffF_invalid];
            c10_poolHffC_invalid ...
                = [c10_poolHffC_invalid; c32_valHffC_invalid; c06_hffC_invalid];
            c10_poolHffCEQ_invalid ...
                = [c10_poolHffCEQ_invalid; c32_valHffCEQ_invalid; c06_hffCEQ_invalid];
        end

        % Surrogate model construction
        c11_surrogateF = trainSurrogate(problem, k, c07_poolX_valid, c08_poolHffF_valid);
        c11_surrogateC = trainSurrogate(problem, k, c07_poolX_valid, c08_poolHffC_valid);
        c11_surrogateCEQ = trainSurrogate(problem, k, c07_poolX_valid, c08_poolHffCEQ_valid);

        % Invalid region model construction
        if size(c09_poolX_invalid,1) == 0
            c12_invalidRegion = [];
        else
            c12_invalidRegion = trainInvalidRegion(problem, c09_poolX_invalid);
        end

        % Generate starting points for optimization using previous results
        if (k == 1), c14_parX = []; end
        c13_startingPoints = compileStartingPoints(problem, ...
            c07_poolX_valid, c08_poolHffF_valid, c08_poolHffC_valid, c08_poolHffCEQ_valid, ...
            c14_parX);

        % Surrogate model-based optimization
        [c14_parX, c15_parSurF, c15_parSurC, c15_parSurCEQ, c16_parSurOut] ...
            = surrogateOptim(problem, c11_surrogateF, c11_surrogateC, c11_surrogateCEQ, ...
                c12_invalidRegion, c13_startingPoints);

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
        c18_parX_valid = c14_parX(idxValid,:);                  % x         valid
        c19_parSurF_valid = c15_parSurF(idxValid,:);            % f(sur)    valid
        c19_parSurC_valid = c15_parSurC(idxValid,:);            % c(sur)    valid
        c19_parSurCEQ_valid = c15_parSurCEQ(idxValid,:);        % ceq(sur)  valid
        c20_parHffF_valid = c17_parHffF(idxValid,:);            % f(hff)    valid
        c20_parHffC_valid = c17_parHffC(idxValid,:);            % c(hff)    valid
        c20_parHffCEQ_valid = c17_parHffCEQ(idxValid,:);        % ceq(hff)  valid
        c21_parX_invalid = c14_parX(idxInvalid,:);              % x         invalid
        c22_parSurF_invalid = c15_parSurF(idxInvalid,:);        % f(sur)    invalid
        c22_parSurC_invalid = c15_parSurC(idxInvalid,:);        % c(sur)    invalid
        c22_parSurCEQ_invalid = c15_parSurCEQ(idxInvalid,:);    % ceq(sur)  invalid
        c23_parHffF_invalid = c17_parHffF(idxInvalid,:);        % f(hff)    invalid
        c23_parHffC_invalid = c17_parHffC(idxInvalid,:);        % c(hff)    invalid
        c23_parHffCEQ_invalid = c17_parHffCEQ(idxInvalid,:);    % ceq(hff)  invalid

        % Validation sampling, high fidelity function evaluations
        [c24_valX, c25_valSurF, c25_valSurC, c25_valSurCEQ] ...
            = samplingValidation(problem, c14_parX, c15_parSurF, c15_parSurC, c15_parSurCEQ);
        [c26_valHffF, c26_valHffC, c26_valHffCEQ] = evaluateHff(problem, c24_valX);
        [idxValid, idxInvalid] ...
            = separateNaN(c24_valX, c26_valHffF, c26_valHffC, c26_valHffCEQ);
        c27_valX_valid = c24_valX(idxValid,:);                  % x         valid
        c28_valSurF_valid = c25_valSurF(idxValid,:);            % f(sur)    valid
        c28_valSurC_valid = c25_valSurC(idxValid,:);            % c(sur)    valid
        c28_valSurCEQ_valid = c25_valSurCEQ(idxValid,:);        % ceq(sur)  valid
        c29_valHffF_valid = c26_valHffF(idxValid,:);            % f(hff)    valid
        c29_valHffC_valid = c26_valHffC(idxValid,:);            % c(hff)    valid
        c29_valHffCEQ_valid = c26_valHffCEQ(idxValid,:);        % ceq(hff)  valid
        c30_valX_invalid = c24_valX(idxInvalid,:);              % x         invalid
        c31_valSurF_invalid = c25_valSurF(idxInvalid,:);        % f(sur)    invalid
        c31_valSurC_invalid = c25_valSurC(idxInvalid,:);        % c(sur)    invalid
        c31_valSurCEQ_invalid = c25_valSurCEQ(idxInvalid,:);    % ceq(sur)  invalid
        c32_valHffF_invalid = c26_valHffF(idxInvalid,:);        % f(hff)    invalid
        c32_valHffC_invalid = c26_valHffC(idxInvalid,:);        % c(hff)    invalid
        c32_valHffCEQ_invalid = c26_valHffCEQ(idxInvalid,:);    % ceq(hff)  invalid

        % Compute convergence and error metrices
        % Eulerian distance (ED) error
        [c33_EDvec, c33_EDavg] = evaluateEulerianDistance(problem, k, ...
                c29_valHffF_valid, c29_valHffC_valid, c29_valHffCEQ_valid, ...
                c28_valSurF_valid, c28_valSurC_valid, c28_valSurCEQ_valid);
        % Hypervolume (HV) and hypervolume ratio (HVR) of predicted Pareto set
        [c34_HVpred, c34_HVRpred] = approxNDHV(problem, c19_parSurF_valid);
        if k == 1
            c34_HVpredHistory = c34_HVpred;
            c34_HVRpredHistory = c34_HVRpred;
        else
            c34_HVpredHistory = [c34_HVpredHistory; c34_HVpred];
            c34_HVRpredHistory = [c34_HVRpredHistory; c34_HVRpred];
        end
        % Normalized residuals of HV and HVR of predicted Pareto set
        if k == 1
            c35_minHVpred = c34_HVpred;
            c35_maxHVpred = c34_HVpred;
            c35_resHVpred = 1;
            c35_resHVpredHistory = 1;
            c35_minHVRpred = c34_HVRpred;
            c35_maxHVRpred = c34_HVRpred;
            c35_resHVRpred = 1;
            c35_resHVRpredHistory = 1;
        else
            c35_minHVpred = min(c35_minHVpred, c34_HVpred);
            c35_maxHVpred = max(c35_maxHVpred, c34_HVpred);
            c35_resHVpred = abs(c34_HVpred - c34_HVpredHistory(end-1,1)) ...
                / abs(c35_maxHVpred - c35_minHVpred);
            c35_resHVpredHistory = [c35_resHVpredHistory; c35_resHVpred];
            c35_minHVRpred = min(c35_minHVRpred, c34_HVRpred);
            c35_maxHVRpred = max(c35_maxHVRpred, c34_HVRpred);
            c35_resHVRpred = abs(c34_HVRpred - c34_HVRpredHistory(end-1,1)) ...
                / abs(c35_maxHVRpred - c35_minHVRpred);
            c35_resHVRpredHistory = [c35_resHVRpredHistory; c35_resHVRpred];
        end
        % Hypervolume (HV) and hypervolume ratio (HVR) of Pareto set of high fidelity results
        HffC = max(c08_poolHffC_valid,[],2);
        HffCEQ = max(abs(c08_poolHffCEQ_valid),[],2);
        iHff = ones(size(c07_poolX_valid,1),1);
        iHff(HffC>problem.control.tolC) = 0;
        iHff(HffCEQ>problem.control.tolCEQ) = 0;
        [~,ndFHF,ndiHF] = ndSort(c07_poolX_valid(iHff,:), c08_poolHffF_valid(iHff,:));
        [c36_HVhff, c36_HVRhff] = approxNDHV(problem, ndFHF(ndiHF == 1,:));
        if k == 1
            c36_HVhffHistory = c36_HVhff;
            c36_HVRhffHistory = c36_HVRhff;
        else
            c36_HVhffHistory = [c36_HVhffHistory; c36_HVhff];
            c36_HVRhffHistory = [c36_HVRhffHistory; c36_HVRhff];
        end
        % Normalized residuals of HV and HVR of Pareto set of high fidelity results
        if k == 1
            c37_minHVhff = c36_HVhff;
            c37_maxHVhff = c36_HVhff;
            c37_resHVhff = 1;
            c37_resHVhffHistory = 1;
            c37_minHVRhff = c36_HVRhff;
            c37_maxHVRhff = c36_HVRhff;
            c37_resHVRhff = 1;
            c37_resHVRhffHistory = 1;
        else
            c37_minHVhff = min(c37_minHVhff, c36_HVhff);
            c37_maxHVhff = max(c37_maxHVhff, c36_HVhff);
            c37_resHVhff = abs(c36_HVhff - c36_HVhffHistory(end-1,1)) ...
                / abs(c37_maxHVhff - c37_minHVhff);
            c37_resHVhffHistory = [c37_resHVhffHistory; c37_resHVhff];
            c37_minHVRhff = min(c37_minHVRhff, c36_HVRhff);
            c37_maxHVRhff = max(c37_maxHVRhff, c36_HVRhff);
            c37_resHVRhff = abs(c36_HVRhff - c36_HVRhffHistory(end-1,1)) ...
                / abs(c37_maxHVRhff - c37_minHVRhff);
            c37_resHVRhffHistory = [c37_resHVRhffHistory; c37_resHVRhff];
        end

        % Compile data to data structure
        scriptCompileDataStructure;
        result.data = [];
        result.data = data;
        result.problem.randomseed = randomseed;

        % Save intermediate result
        if (problem.nested.outeriter == 0) && (problem.nested.innercase == 0)
            save(fullfile(problem.control.solpath, ...
                [problem.control.case, '_iter', num2str(k,'%04d'), '.mat']), ...
                'result', '-v7.3');
        else
            save(fullfile(problem.control.solpath, ...
                [problem.control.case, ...
                '_outeriter', num2str(problem.nested.outeriter,'%04d'), ...
                '_innercase', num2str(problem.nested.innercase,'%04d'), ...
                '_iter', num2str(k,'%04d'), '.mat']), ...
                'result', '-v7.3');
        end

        % Stopping criteria evaluation
        if evaluateStopCriteria(), break; end
    end

end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
