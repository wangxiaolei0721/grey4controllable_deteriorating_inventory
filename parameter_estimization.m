%% clear data and figure
clc;
clear;
close all;
%% load data
load(".\data\external_settings.mat")
load(".\data\simulated_data.mat")
%% data partitioning
train_length = 0.6 * m;
% simulated time 1,2,3 for parameter estimation
time_train=time_simu(1:train_length);
demand_train=demand_simu(1:train_length);
level_train=level_simu(1:train_length);
xi_vector_train = xi_vector(1:train_length);
p_vector_train = p_vector(1:train_length);
Q_vector_train = Q_vector(1:train_length);
%% parameter estimation
% the initial value of theta
inventory_pars0 = inventory_pars_initial(xi_vector_train,demand_train,Q_vector_train,level_train);
% inventory_pars0=[0.1;4];
inventory_var0 = inventory_pars2var(xi_vector_train,demand_train,Q_vector_train,level_train,inventory_pars0);
% the initial demand parameter correponding to deterioration parameter
[demand_pars0,demand_var0] = inventory_pars2demand_pars(time0,time_train,xi_vector_train,p_vector_train,demand_train,inventory_pars0);
% weight definition
weight_initial=[1/demand_var0;1/inventory_var0];
% maximum number of iterations
max_iter=1000;
% Iteratively Reweighed Least Squares algorithm for parameter estimation
% previous deterioration parameter
inventory_pars_prev = inventory_pars0;
% record the estimated parameter values for each iteration
history = zeros(2,max_iter); 
history(:,1) = inventory_pars0;
weight=weight_initial;
% lsqnonlin function setting
opt_options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',10,'FunctionTolerance',1e-5,'StepTolerance',1e-4);
for iter = 2:max_iter
    disp(iter)
    % residual function of lambda under multiple orders 
    minobjfun = @(inventory_pars) inventory_pars_res(time0,time_train,xi_vector_train,p_vector_train,Q_vector_train,demand_train,level_train,weight,inventory_pars);
    % next theta
    inventory_pars_next = lsqnonlin(minobjfun,inventory_pars_prev,0,0.2,opt_options);
    % the residual variance of inventory regression equation based on theta
    inventory_var = inventory_pars2var(xi_vector_train,demand_train,Q_vector_train,level_train,inventory_pars_next);
    % the initial values of alpha and beta based on theta
    [~,demand_var] = inventory_pars2demand_pars(time0,time_train,xi_vector_train,p_vector_train,demand_train,inventory_pars_next);
    % weight
    weight=[1/demand_var;1/inventory_var];
    % record parameter estimates for the current iteration
    history(:,iter) = inventory_pars_next;
    inventory_pars_prev=inventory_pars_next;
end
inventory_pars=inventory_pars_next;
%  the estimated values of alpha and beta based on theta
[demand_pars,demand_var_train] = inventory_pars2demand_pars(time0,time_train,xi_vector_train,p_vector_train,demand_train,inventory_pars);
inventory_var_train = inventory_pars2var(xi_vector_train,demand_train,Q_vector_train,level_train,inventory_pars);
alpha_estimate=demand_pars(1);
beta_estimate=demand_pars(2);
lambda0_estimate=inventory_pars(1);
rho_estimate=inventory_pars(2);
disp(inventory_pars - inventory_pars0)
save(".\data\parameter.mat","lambda0_estimate","rho_estimate","alpha_estimate","beta_estimate")
%% fitted level
% test
time_test=time_simu((train_length+1):end);
demand_test=demand_simu((train_length+1):end);
level_test=level_simu((train_length+1):end);
xi_vector_test = xi_vector((train_length+1):end);
p_vector_test = p_vector((train_length+1):end);
Q_vector_test = Q_vector((train_length+1):end);
% test error
demand_var_test = demand_pars2demand_vars(time0,time_test,xi_vector_test,p_vector_test,demand_test,inventory_pars,demand_pars);
inventory_var_test = inventory_pars2var(xi_vector_test,demand_test,Q_vector_test,level_test,inventory_pars);
% plot
time_fit = {};
demand_fit = {};
level_fit = {};
for i = 1:m
    [time_fit_i,demand_fit_i,level_fit_i] = inventory_level(alpha_estimate,beta_estimate,p_vector(i),lambda0_estimate,rho_estimate,xi_vector(i),time0,delta_t,Q_vector(i));
    time_fit{i} = time_fit_i; % [time0;time_fit_i];
    demand_fit{i} = demand_fit_i;
    level_fit{i} = level_fit_i; % [Q_vector(i);level_fit_i];
end
%% plot
% load data
load(".\data\true_data.mat")
% demand plot
fdemand=figure('unit','centimeters','position',[5,5,40,20],'PaperPosition',[5,5,40,20],'PaperSize',[40,20]);
tiledlayout(2,m/2,'Padding','Compact');
% plot time vs level
finvertory=figure('unit','centimeters','position',[5,5,40,20],'PaperPosition',[5,5,40,20],'PaperSize',[40,20]);
tiledlayout(2,m/2,'Padding','Compact');
for i = 1:m
    figure(fdemand)
    nexttile
    plot(time_true{i},demand_true{i},'LineWidth',1.5,'Color',[0 0.4470 0.7410],'Marker','^','MarkerSize',8)
    hold on
    plot(time_simu{i},demand_simu{i},'LineWidth',1.5,'LineStyle','--','Color',[0.9290 0.6940 0.1250],'Marker','o','MarkerSize',8)
    plot(time_fit{i},demand_fit{i},'LineWidth',1.5,'Color',[0.8500 0.3250 0.0980],'Marker','square','MarkerSize',8)
    xlabel({'Day'},'FontSize',12)
    xlim([min(time_simu{i}),max(time_simu{i})])
    ylabel(['Demand'],'FontSize',12)
    title(strcat("(",char(96 + i),") The ", num2str(i),"th ordering cycle"),'FontSize',14)
    set(gca,'FontName','Book Antiqua','FontSize',10)
    if i==10
        legend(["Actual demand","Simulated demand","Fitted demand"],'location','northeast','FontSize',10,'NumColumns',1)
    end
    figure(finvertory)
    nexttile
    plot(time_true{i},level_true{i},'LineStyle','-','LineWidth',1.5,'Color',[0 0.4470 0.7410],'Marker','^','MarkerSize',6)
    hold on
    plot(time_simu{i},level_simu{i},'LineStyle','--','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250],'Marker','o','MarkerSize',6)
    plot(time_fit{i},level_fit{i},'LineStyle','-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980],'Marker','square','MarkerSize',6)
    xlabel({'Day'},'FontSize',12)
    xlim([min(time_simu{i}),max(time_simu{i})])
    ylabel(['Inventory level'],'FontSize',12)
    title(strcat("(",char(96 + i),") The ", num2str(i),"th ordering cycle"),'FontSize',14)
    set(gca,'FontName','Book Antiqua','FontSize',10)
    if i==10
        legend(["Actual inventory level","Simulated inventory level","Fitted inventory level"],'location','northeast','FontSize',8,'NumColumns',1) % ,"Fitted inventory level"
    end
end
% save figure
savefig(fdemand,'.\figure\simulation_demand.fig')
exportgraphics(fdemand,'.\figure\simulation_demand.pdf')
savefig(finvertory,'.\figure\simulation_level.fig')
exportgraphics(finvertory,'.\figure\simulation_level.pdf')


