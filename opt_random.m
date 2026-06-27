% clear data and figure
clc;
clear;
close all;
tic
%% model setting
% load estimated parameters
load(".\data\parameter.mat")
%% economic order quantity
c=8.0;
h=0.2;
K=100;
p=12;
%% profit different xi for profit
xi_interval=[0 1];
% cycle interval
T_interval=[1 7];
profit_fd = @(T,xi) profit(K,c,h,p,time0,pars_true,T,xi);
figure('unit','centimeters','position',[25,5,15,10],'PaperPosition',[25,5,15,10],'PaperSize',[15,10]);
fsurf(profit_fd,[T_interval,xi_interval])
hold on;
%% solve
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
T_opt  = x_opt(1);
xi_opt = x_opt(2);
profit_opt = -fval;
% T_opt
Q_opt=pars_true(4)*T_opt;
plot3(T_opt,xi_opt,profit_opt,'Marker','hexagram','MarkerSize',12,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0.2 0.2],'LineStyle','none')
xlabel(['Order cycle'],'FontSize',12)
ylabel(['Investment'],'FontSize',12)
zlabel(['Profit'],'FontSize',12)
% ylim([-70,45])
set(gca,'FontName','Book Antiqua','FontSize',12)
title(['Average profit'],'FontSize',10)
%% multi pars
T_opt_ini = [];
xi_opt_ini = [];
Q_opt_ini = [];
profit_opt_ini = [];
%
T_opt_est = [];
xi_opt_est = [];
Q_opt_est = [];
profit_opt_est = [];
for i =1:1000
    %
    pars=pars_ini(:,i);
    %% solve
    obj = @(x) -profit(K,c,h,p,time0,pars,x(1),x(2));
    problem = createOptimProblem('fmincon', ...
        'objective', obj, ...
        'x0', [3 0.2], ...
        'lb', [1 0], ...
        'ub', [7 0.5]);
    ms = MultiStart;
    [x_opt,fval] = run(ms,problem,20);
    T_opt_ini(i,1)  = x_opt(1);
    xi_opt_ini(i,1) = x_opt(2);
    profit_opt_ini(i,1) = -fval;
    % T_opt
    Q_opt_ini(i,1) = pars(4)*T_opt_ini(i);
    %
    pars=pars_est(:,i);
    %% solve
    obj = @(x) -profit(K,c,h,p,time0,pars,x(1),x(2));
    problem = createOptimProblem('fmincon', ...
        'objective', obj, ...
        'x0', [3 0.2], ...
        'lb', [1 0], ...
        'ub', [7 0.5]);
    ms = MultiStart;
    [x_opt,fval] = run(ms,problem,20);
    T_opt_est(i,1)  = x_opt(1);
    xi_opt_est(i,1) = x_opt(2);
    profit_opt_est(i,1) = -fval;
    % T_opt
    Q_opt_est(i,1) =pars(4)*T_opt_est(i);
end
% boxplot
labels = ["initial ests", "Final ests"];
fig_opt = figure('unit','centimeters','position',[10,10,30,10],'PaperPosition',[0, 0, 30,10],'PaperSize',[30,10]);
ax = subplot(1,2,1);
order_compare = [T_opt_ini, T_opt_est];
bp = boxplot(order_compare,'Labels',labels);
hold on
yline(T_opt,'r--','LineWidth',1.5);
set(findobj(gca,'Type','Line'),'LineWidth',1.5);
set(gca,'FontName','Book Antiqua','FontSize',12)
title(ax, "(a) Replenishment cycle", 'FontSize',14);
ax = subplot(1,2,2);
invest_compare = [xi_opt_ini, xi_opt_est];
bp = boxplot(invest_compare,'Labels',labels);
hold on
yline(xi_opt,'r--','LineWidth',1.5);
set(findobj(gca,'Type','Line'),'LineWidth',1.5);
set(gca,'FontName','Book Antiqua','FontSize',12)
title(ax, "(b) Preservation investment", 'FontSize',14);
% save(".\data\parameter.mat","lambda_estimate","eta_estimate","rho_estimate","d_estimate")
savefig(fig_opt,'.\figure\simulation_opt.fig')
exportgraphics(fig_opt,'.\figure\simulation_opt.pdf')
toc





