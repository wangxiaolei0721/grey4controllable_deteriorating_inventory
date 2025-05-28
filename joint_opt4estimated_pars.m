% clear data and figure
clc;
clear;
close all;
%% model setting
% load estimated parameters
load(".\data\parameter.mat")
%% economic order quantity
c=3.0;
h=0.1;
K=50;
% price interval based on estimates
p_fit_interval=[c alpha_estimate/beta_estimate];
xi_interval=[0 2];
% cycle interval
T_interval=[1 7];
%% profit different xi for true profit
figure('unit','centimeters','position',[5,5,30,15],'PaperPosition',[5,5,30,15],'PaperSize',[30,15]);
tiledlayout(1,2,'Padding','Compact');
nexttile
xi1=0;
xi2=0.2;
lambda1 = lambda0_estimate*exp(-rho_estimate*xi1);
lambda2 = lambda0_estimate*exp(-rho_estimate*xi2);
profit_appro_fd0 = @(p,T) profit_appro(alpha_estimate,beta_estimate,p,lambda0_estimate,rho_estimate,xi1,c,h,K,T);
fsurf(profit_appro_fd0,[p_fit_interval,T_interval])
xlabel({'Price'},'FontSize',12)
ylabel(['Order cycle/day'],'FontSize',12)
zlabel(['Profit/unit currency'],'FontSize',12)
zlim([-70,45])
set(gca,'FontName','Book Antiqua','FontSize',12)
title(['(a) Approximate profit surface without considering preservation investment'],'FontSize',10)
nexttile
profit_appro_fd1 = @(p,T) profit_appro(alpha_estimate,beta_estimate,p,lambda0_estimate,rho_estimate,xi2,c,h,K,T);
fsurf(profit_appro_fd1,[p_fit_interval,T_interval])
xlabel({'Price'},'FontSize',12)
ylabel(['Order cycle/day'],'FontSize',12)
zlabel(['Profit/unit currency'],'FontSize',12)
zlim([-70,45])
set(gca,'FontName','Book Antiqua','FontSize',12)
title(['(b) Approximate profit surface when preservation investment \xi_2=0.2'],'FontSize',10)
savefig(gcf,'.\figure\fit_analysis_xi.fig')
exportgraphics(gcf,'.\figure\fit_analysis_xi.pdf')
% together
figure('unit','centimeters','position',[5,5,15,15],'PaperPosition',[5,5,15,15],'PaperSize',[15,15]);
fsurf(profit_appro_fd0,[p_fit_interval,T_interval])
hold on
fsurf(profit_appro_fd1,[p_fit_interval,T_interval])
xlabel({'Price'},'FontSize',12)
ylabel(['Order cycle/day'],'FontSize',12)
zlabel(['Profit/unit currency'],'FontSize',12)
zlim([-70,45])
legend(["Preservation investment \xi_1=0","Preservation investment \xi_2=0.2"],'location','northeast','FontSize',8,'NumColumns',1)
set(gca,'FontName','Book Antiqua','FontSize',12)
savefig(gcf,'.\figure\fit_analysis_xi_together.fig')
exportgraphics(gcf,'.\figure\fit_analysis_xi_together.pdf')
%% solve
syms p xi;
% T_opt
lambda=lambda0_estimate*exp(-rho_estimate*xi);
denominator=(alpha_estimate-beta_estimate*p)*(p*lambda+h+xi);
T_p_xi=sqrt(2*K/denominator);
profit_syms = profit_appro(alpha_estimate,beta_estimate,p,lambda0_estimate,rho_estimate,xi,c,h,K,T_p_xi);
profit_der_p=diff(profit_syms,p);
profit_der_xi=diff(profit_syms,xi);
eq1 = profit_der_p == 0;
eq2 = profit_der_xi == 0;
sol = vpasolve([eq1, eq2], [p, xi],[p_fit_interval;xi_interval]);
p_opt  = double(sol.p);
xi_opt  = double(sol.xi);
% T_opt
lambda=lambda0_estimate*exp(-rho_estimate*xi_opt);
denominator=(alpha_estimate-beta_estimate*p_opt)*(p_opt*lambda+h+xi_opt);
T_opt=sqrt(2*K/denominator);
Q_opt=(alpha_estimate-beta_estimate*p_opt)*T_opt;
profit_opt = profit_appro(alpha_estimate,beta_estimate,p_opt,lambda0_estimate,rho_estimate,xi_opt,c,h,K,T_opt);
% plot
figure('unit','centimeters','position',[5,5,15,15],'PaperPosition',[5,5,15,15],'PaperSize',[15,15]);
profit_fd = @(p,xi) profit_appro(alpha_estimate,beta_estimate,p,lambda0_estimate,rho_estimate,xi,c,h,K,T_opt);
fsurf(profit_fd,[p_fit_interval,xi_interval])
hold on
plot3(p_opt,xi_opt,profit_opt,'Marker','hexagram','MarkerSize',12,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0.2 0.2],'LineStyle','none')
set(gca,'FontName','Book Antiqua','FontSize',12)
xlabel({'Price'},'FontSize',12)
ylabel(['Investment'],'FontSize',12)
zlabel(['Profit/unit currency'],'FontSize',12)

% save figure
savefig(gcf,'.\figure\fit_opt.fig')
exportgraphics(gcf,'.\figure\fit_opt.pdf')




