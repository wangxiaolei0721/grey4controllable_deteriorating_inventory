% clear data and figure
clc;
clear;
close all;
%% model setting
% load estimated parameters
load(".\data\parameter_case.mat")
%% economic order quantity
c=8.0;
h=0.2;
K=100;
p=12;
% profit
xi_interval=[0 1];
% cycle interval
T_interval=[1 7];
figure('unit','centimeters','position',[5,5,30,15],'PaperPosition',[5,5,30,15],'PaperSize',[30,15]);
tiledlayout(1,2,'Padding','Compact');
%%
nexttile
profit_fd = @(T,xi) profit(K,c,h,p,time0,pars_true,T,xi);
fsurf(profit_fd,[T_interval,xi_interval])
hold on;
% solve
obj = @(x) -profit(K,c,h,p,time0,pars_true,x(1),x(2));
% x0 = [3, 0.2];
% lb = [1, 0];
% ub = [7, 0.5];
% options = optimoptions('fmincon','Display','iter');
% [x_opt, fval] = fmincon(obj,x0,[],[],[],[],lb,ub,[],options);
% multi initial point
problem = createOptimProblem('fmincon', ...
    'objective', obj, ...
    'x0', [3 0.2], ...
    'lb', [1 0], ...
    'ub', [7 0.5]);
ms = MultiStart;
[x_opt,fval] = run(ms,problem,20);
T_opt_true  = x_opt(1);
xi_opt_true = x_opt(2);
profit_opt_true = -fval;
% T_opt
Q_opt_true = pars_true(4)*T_opt_true;
plot3(T_opt_true,xi_opt_true,profit_opt_true,'Marker','hexagram','MarkerSize',12,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0.2 0.2],'LineStyle','none')
xlabel(['Replenishment cycle'],'FontSize',12)
ylabel(['Preservation investment'],'FontSize',12)
zlabel(['Profit'],'FontSize',12)
% ylim([-70,45])
set(gca,'FontName','Book Antiqua','FontSize',12)
title(['(a) Profit for true parameters'],'FontSize',14)
%%
nexttile
profit_fd = @(T,xi) profit(K,c,h,p,time0,pars_est,T,xi);
fsurf(profit_fd,[T_interval,xi_interval])
hold on;
% solve
obj = @(x) -profit(K,c,h,p,time0,pars_est,x(1),x(2));
% x0 = [3, 0.2];
% lb = [1, 0];
% ub = [7, 0.5];
% options = optimoptions('fmincon','Display','iter');
% [x_opt, fval] = fmincon(obj,x0,[],[],[],[],lb,ub,[],options);
% multi 初始值
problem = createOptimProblem('fmincon', ...
    'objective', obj, ...
    'x0', [3 0.2], ...
    'lb', [1 0], ...
    'ub', [7 0.5]);
ms = MultiStart;
[x_opt,fval] = run(ms,problem,20);
T_opt_est  = x_opt(1);
xi_opt_est = x_opt(2);
profit_opt_est = -fval;
% T_opt
Q_opt_est = pars_est(4)*T_opt_est;
plot3(T_opt_est,xi_opt_est,profit_opt_est,'Marker','hexagram','MarkerSize',12,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0.2 0.2],'LineStyle','none')
xlabel(['Replenishment cycle'],'FontSize',12)
ylabel(['Preservation investment'],'FontSize',12)
zlabel(['Profit'],'FontSize',12)
% ylim([-70,45])
set(gca,'FontName','Book Antiqua','FontSize',12)
title(['(b) Profit for estimated parameters'],'FontSize',14)
% save figure
savefig(gcf,'.\figure\simulation_opt_case.fig')
exportgraphics(gcf,'.\figure\simulation_opt_case.pdf')


