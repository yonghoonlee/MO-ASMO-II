set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultColorbarTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultPolaraxesTickLabelInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultTextarrowshapeInterpreter','latex');
set(0,'DefaultTextboxshapeInterpreter','latex');
close all;
load('solution/Lee2017a_NonlconCheapresults.mat');
fg = figure('Color',[1 1 1]);
fg.Position = [100 100 560 420];
data = resultNSGA2.DO.fopt;
ph = plot(data(:,1),data(:,2),'ko','MarkerFaceColor',[0 0 0.9]);
ax = gca;
xlabel('$f_1(\mathbf{x})$');
ylabel('$f_2(\mathbf{x})$');
grid on;
ax.GridLineStyle = ':';
ax.GridColor = [0 0 0];
ax.GridAlpha = 1;
ax.LineWidth = 1.5;
lh = legend('Pareto front');
ax.FontSize = 13;
lh.FontSize = 13;
lh.Interpreter = 'tex';
export_fig 'Lee2017a' -pdf

%---------

[X,Y] = meshgrid(-5:0.1:5, -5:0.1:5);
XY = [reshape(X,numel(X),1) reshape(Y,numel(Y), 1)];
F1 = (3*sin(2.5*XY(:,1)) - 2*XY(:,1)).*cos(XY(:,2)).*exp(-(1e-3)*XY(:,2).^2);
F2 = (3/8)*(0.3*abs(XY(:,1))).^(19/25).*(sin(5*XY(:,2)).*XY(:,2));
FV = [F1, F2];
FM1 = reshape(F1, size(X,1), size(X,2));
FM2 = reshape(F2, size(X,1), size(X,2));
CV = (100*(XY(:,2) - XY(:,1).^2).^2+(XY(:,1) - 1).^2) - 1;
CM = reshape(CV, size(X,1), size(X,2));
fg1 = figure('Color',[1 1 1]);
fg1.Position = [400 120 560 420];
contourf(X,Y,FM1);
fg2 = figure('Color',[1 1 1]);
fg2.Position = [700 140 560 420];
contourf(X,Y,FM2);
fg3 = figure('Color',[1 1 1]);
fg3.Position = [1000 160 560 420];
contourf(X,Y,CM);

%---------

% Assume nonlinear constraint is cheap
load('solution/Lee2017a_NonlconCheapresults.mat');
dMOC = resultMOASMO.data;
dGAC = resultNSGA2.DO;

% Assume nonlinear constraint is expensive
load('solution/Lee2017a_NonlconExpresults.mat');
dMOE = resultMOASMO.data;
dGAE = resultNSGA2.DO;

dataC = dMOC.c15_parSurF{1,1};
dataE = dMOE.c15_parSurF{1,1};

fg4 = figure('Color',[1 1 1]);
fg4.Position = [1300 180 560 420];
plot(dataC(:,1), dataC(:,2), 'ko', 'MarkerFaceColor', 'r'); hold on;
plot(dataE(:,1), dataE(:,2), 'kx'); hold on;
