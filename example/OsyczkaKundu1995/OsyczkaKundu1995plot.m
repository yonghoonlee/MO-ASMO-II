[path1,file1] = fileparts(mfilename('fullpath')); [path2,file2] = fileparts(path1);
path3 = fullfile(path1,'solution'); file3 = [file2,'results.mat'];
load(fullfile(path3,file3));

xopt1 = resultMOASMO.data.c07_poolX_valid{1,1};
fopt1 = resultMOASMO.data.c08_poolHffF_valid{1,1};
copt1 = resultMOASMO.data.c08_poolHffC_valid{1,1};
ceqopt1 = resultMOASMO.data.c08_poolHffCEQ_valid{1,1};
if (size(copt1,2) ~= 0)
    copt1 = max(copt1, [], 2);
else
    copt1 = zeros(size(fopt1,1),1);
end
if (size(ceqopt1,2) ~= 0)
    ceqopt1 = max(ceqopt1.^2, [], 2);
else
    ceqopt1 = zeros(size(fopt1,1),1);
end

ioptm = (((copt1 - resultNSGA2.problem.control.tolC) <= 0) ...
    & ((sqrt(ceqopt1) - resultNSGA2.problem.control.tolCEQ) <= 0));
ioptm = enforceIndexLincon(problem, ioptm, xopt1);
xopt1 = xopt1(ioptm, :);
fopt1 = fopt1(ioptm, :);
[xopt1, fopt1, iopt1] = ndSort(xopt1, fopt1);

cm = plasma(max(iopt1)+5);
fg = figure();
for idx = 1:max(iopt1)
    plot(fopt1(iopt1==idx,1),fopt1(iopt1==idx,2),'o','Color',cm(idx+5,:)); hold on;
end
xopt2 = resultNSGA2.DO.xopt;
fopt2 = resultNSGA2.DO.fopt;
[xopt2, fopt2, iopt2] = ndSort(xopt2, fopt2);
plot(fopt2(iopt2==1,1),fopt2(iopt2==1,2),'.','Color',cm(1,:)); hold on;