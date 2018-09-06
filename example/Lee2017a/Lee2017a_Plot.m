close all;
clear; clc;

set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultColorbarTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultPolaraxesTickLabelInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultTextarrowshapeInterpreter','latex');
set(0,'DefaultTextboxshapeInterpreter','latex');
set(0,'DefaultAxesFontSize', 14);
set(0,'DefaultColorbarFontSize', 12);
set(0,'DefaultLegendFontSize', 12);
set(0,'DefaultPolaraxesFontSize', 14);
set(0,'DefaultTextFontSize', 14);
set(0,'DefaultTextarrowshapeFontSize', 14);
set(0,'DefaultTextboxshapeFontSize', 14);
set(0,'DefaultUibuttongroupFontSize', 14);
set(0,'DefaultUicontrolFontSize', 14);
set(0,'DefaultUipanelFontSize', 14);
set(0,'DefaultUitableFontSize', 14);
set(0,'defaultFigurePosition', [680 620 560 320]);

resultCheapNonlcon = load('/Users/yonghoonlee/Dropbox/GitRepositories/MO-ASMO-II/example/Lee2017a/solution/Lee2017a_CheapNonlcon_ECresults.mat');
resultExpNonlcon = load('/Users/yonghoonlee/Dropbox/GitRepositories/MO-ASMO-II/example/Lee2017a/solution/Lee2017a_ExpNonlcon_ECresults.mat');

figure('Color',[1 1 1],'Position',[480 390 680 500]);
subplot(2,2,1);
ph1 = plot(resultCheapNonlcon.resultMOASMO.fopt(:,1), resultCheapNonlcon.resultMOASMO.fopt(:,2), ...
    's-','Color',[1 0.25 0.25],'MarkerFaceColor',[1 0.25 0.25],'MarkerEdgeColor','none','MarkerSize',10); hold on;
ph2 = plot(resultExpNonlcon.resultMOASMO.fopt(:,1), resultExpNonlcon.resultMOASMO.fopt(:,2), ...
    'd-','Color',[0 0 0.75],'MarkerFaceColor',[0 0 0.75],'MarkerEdgeColor','none','MarkerSize',8); hold on;
lgd1 = legend([ph1 ph2], {'constrained sampling','standard sampling'}, ...
    'FontSize', 20, 'Box', 'off');
thx1 = xlabel('$f_1$'); thy1 = ylabel('$f_2$');
tht1 = title('(a)'); tht1.Position = [-0.5165 -0.9029];
ax = gca; ax.FontSize = 20;
ax.XLim = [-1 7]; ax.XTick = -1:7; ax.YLim = [-0.9 0.9]; ax.YTick = -0.9:0.3:0.3;
lgd1.Position = [0.12835 0.824 0.33195 0.105];

subplot(2,2,2);
xpoolall = [resultCheapNonlcon.resultMOASMO.data.c07_poolX_valid{1,1};
            resultCheapNonlcon.resultMOASMO.data.c27_valX_valid{1,1}];
fpoolall = [resultCheapNonlcon.resultMOASMO.data.c08_poolHffF_valid{1,1};
            resultCheapNonlcon.resultMOASMO.data.c29_valHffF_valid{1,1}];
cpoolall = [resultCheapNonlcon.resultMOASMO.data.c08_poolHffC_valid{1,1};
            resultCheapNonlcon.resultMOASMO.data.c29_valHffC_valid{1,1}];
ceqpoolall = [resultCheapNonlcon.resultMOASMO.data.c08_poolHffCEQ_valid{1,1};
              resultCheapNonlcon.resultMOASMO.data.c29_valHffCEQ_valid{1,1}];
ipoolall = ones(size(fpoolall, 1), 1);
if size(cpoolall, 2) > 0
    cpoolall = max(cpoolall, [], 2);
    icpool = zeros(size(ipoolall));
    icpool((cpoolall - resultCheapNonlcon.resultMOASMO.problem.control.tolC) <= 0) = 1;
    ipoolall = ipoolall & icpool;
end
if size(ceqpoolall, 2) > 0
    ceqpoolall = max(sqrt(ceqpoolall.^2), [], 2);
    iceqpool = zeros(size(ipoolall));
    iceqpool((ceqpoolall - resultCheapNonlcon.resultMOASMO.problem.control.tolCEQ) <= 0) = 1;
    ipoolall = ipoolall & iceqpool;
end
ipoolall = logical(ipoolall);
xfea = xpoolall(ipoolall, :);
ffea = fpoolall(ipoolall, :);
xinfea = xpoolall(~ipoolall, :);
finfea = fpoolall(~ipoolall, :);
ph3 = plot(finfea(:,1), finfea(:,2), ...
    '.','MarkerFaceColor','none','Color',[1 0.5 0.5],'MarkerSize',10); hold on;
ph4 = plot(ffea(:,1), ffea(:,2), ...
    's','MarkerFaceColor',[1 0.15 0.15],'MarkerEdgeColor','none','MarkerSize',10); hold on;

xpoolall = [resultExpNonlcon.resultMOASMO.data.c07_poolX_valid{1,1};
            resultExpNonlcon.resultMOASMO.data.c27_valX_valid{1,1}];
fpoolall = [resultExpNonlcon.resultMOASMO.data.c08_poolHffF_valid{1,1};
            resultExpNonlcon.resultMOASMO.data.c29_valHffF_valid{1,1}];
cpoolall = [resultExpNonlcon.resultMOASMO.data.c08_poolHffC_valid{1,1};
            resultExpNonlcon.resultMOASMO.data.c29_valHffC_valid{1,1}];
ceqpoolall = [resultExpNonlcon.resultMOASMO.data.c08_poolHffCEQ_valid{1,1};
              resultExpNonlcon.resultMOASMO.data.c29_valHffCEQ_valid{1,1}];
ipoolall = ones(size(fpoolall, 1), 1);
if size(cpoolall, 2) > 0
    cpoolall = max(cpoolall, [], 2);
    icpool = zeros(size(ipoolall));
    icpool((cpoolall - resultExpNonlcon.resultMOASMO.problem.control.tolC) <= 0) = 1;
    ipoolall = ipoolall & icpool;
end
if size(ceqpoolall, 2) > 0
    ceqpoolall = max(sqrt(ceqpoolall.^2), [], 2);
    iceqpool = zeros(size(ipoolall));
    iceqpool((ceqpoolall - resultExpNonlcon.resultMOASMO.problem.control.tolCEQ) <= 0) = 1;
    ipoolall = ipoolall & iceqpool;
end
ipoolall = logical(ipoolall);
xfea = xpoolall(ipoolall, :);
ffea = fpoolall(ipoolall, :);
xinfea = xpoolall(~ipoolall, :);
finfea = fpoolall(~ipoolall, :);
ph5 = plot(finfea(:,1), finfea(:,2), ...
    '.','MarkerFaceColor','none','Color',[0.15 0.15 1],'MarkerSize',10); hold on;
ph6 = plot(ffea(:,1), ffea(:,2), ...
    'd','MarkerFaceColor',[0 0 0.75],'MarkerEdgeColor','none','MarkerSize',10); hold on;

thx2 = xlabel('$f_1$'); thy2 = ylabel('$f_2$');
ax = gca; ax.FontSize = 20;
ax.XLim = [-12 12]; ax.XTick = -12:4:12; ax.YLim = [-2 3]; ax.YTick = -2:4;
tht2 = title('(b)'); tht2.Position = [10.4900 -2.0126];
lgd2 = legend([ph6, ph5], {'Feasible point','Infeasible point'},'FontSize',20,'box','off');
lgd2.Position = [0.5560 0.8260 0.2698 0.1050];


% save
if resultExpNonlcon.resultMOASMO.problem.control.plotexport
    if (exist(resultExpNonlcon.resultMOASMO.problem.control.plotpath) == 0)
        mkdir(resultExpNonlcon.resultMOASMO.problem.control.plotpath)
    end
    fnm = fullfile(resultExpNonlcon.resultMOASMO.problem.control.plotpath, ...
        'Lee2017a_COMPARE_plot_COMP_CONV');
    eval(['export_fig ''', fnm, ''' -pdf']);
end