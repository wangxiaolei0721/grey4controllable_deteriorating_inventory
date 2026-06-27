% clear data and figure
clc;
clear;
close all;
%% model setting
% deterioration parameters
lambda = 0.05;
eta = 2;
%
rho=2;
% demand parameters
demand = 50;
pars_true = [lambda;eta;rho;demand];
%% simulation parameters
% the time of order arrival
time0=0;
% the time resolution
delta_t=1;
% standard deviation of error in demand regression equation54
std_dev=2;
% order cycles
m=10;
% data storage
pars_ini=[];
pars_est=[];
%
demand_error_ini = [];
level_diff_error_ini = [];
error_ini = [];
%
demand_error_est = [];
level_diff_error_est = [];
level_error_est = [];
%
for j=1:1000
    % random seed
    rng(j)
    % the order quantity
    Q_vector=300+randi([0,60],m,1);
    % the sales price
    % the preservation cost
    xi_vector = 0.2*rand(m,1); % 0.5*rand(m,1)
    %% initialization of data storage
    time_true = {};
    time_true_t0 = {};
    demand_true = {};
    level_diff_true = {};
    level_true = {};
    level_true_Q = {};
    %
    time_simu = {};
    time_simu_t0 = {};
    demand_simu = {};
    level_simu_Q = {};
    % generate the inventory levels
    for i = 1:m
        % true level
        [time_true_i,demand_true_i,level_diff_true_i,level_true_i] = inventory_level(demand,lambda,eta,rho,xi_vector(i),time0,delta_t,Q_vector(i));
        time_true{i}=time_true_i;
        time_true_t0{i} = [time0;time_true_i];
        demand_true{i} = demand_true_i;
        level_diff_true{i}=level_diff_true_i;
        level_true{i}=level_true_i;
        level_true_Q{i}=[Q_vector(i);level_true_i];
        % simulated level
        [time_simu_i,demand_simu_i,level_diff_simu_i,level_simu_i] = inventory_level_simulation(demand,std_dev,lambda,eta,rho,xi_vector(i),time0,delta_t,Q_vector(i));
        time_simu{i} = time_simu_i;
        time_simu_t0{i} = [time0;time_simu_i];
        demand_simu{i} = demand_simu_i;
        level_diff_simu{i}=level_diff_simu_i;
    end
    %% cumulative generation
    for i = 1:m
        time_i=time_simu_t0{i};
        time_i_diff=diff(time_i);
        level_diff_i=level_diff_simu{i};
        level_simu_i = [Q_vector(i);Q_vector(i) + cumsum(level_diff_i.*time_i_diff)];
        level_simu_Q{i}=level_simu_i;
    end
    %% data partitioning
    train_length = 0.6 * m;
    % simulated time 1,2,3 for parameter estimation
    time_train=time_simu(1:train_length);
    time_t0_train=time_simu_t0(1:train_length);
    demand_train=demand_simu(1:train_length);
    level_diff_train=level_diff_simu(1:train_length);
    level_train=level_simu_Q(1:train_length);
    xi_vector_train = xi_vector(1:train_length);
    %% parameter estimation
    % the initial value of theta time0,time_train,xi_vector_train,demand_train,level_diff_train,level_train
    pars0 = inventory_pars_initial(time_t0_train,xi_vector_train,demand_train,level_diff_train,level_train);
    inventory_var0 = inventory_pars2var(time_t0_train,xi_vector_train,demand_train,level_diff_train,level_train,pars0);
    % the initial demand parameter correponding to deterioration parameter
    [demand_par0,demand_var0] = inventory_pars2demand_pars(time_t0_train,xi_vector_train,demand_train,pars0);
    % weight definition
    weight_initial=[1/demand_var0;1/inventory_var0];
    % maximum number of iterations
    max_iter=100;
    % Iteratively Reweighed Least Squares algorithm for parameter estimation
    % previous deterioration parameter
    pars_prev = pars0;
    % record the estimated parameter values for each iteration
    history = zeros(3,max_iter);
    history(:,1) = pars0;
    weight=weight_initial;
    % lsqnonlin function setting
    opt_options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',10,'FunctionTolerance',1e-5,'StepTolerance',1e-4);
    for iter = 2:max_iter
        %         disp(iter)
        % residual function of lambda under multiple orders
        minobjfun = @(inventory_pars) inventory_pars_res(time_t0_train,xi_vector_train,demand_train,level_diff_train,level_train,weight,inventory_pars);
        % next theta
        pars_next = lsqnonlin(minobjfun,pars_prev,[0,0,0],[0.5,3,10], opt_options);
        % the residual variance of inventory regression equation based on theta
        inventory_var = inventory_pars2var(time_t0_train,xi_vector_train,demand_train,level_diff_train,level_train,pars_next);
        % the initial values of alpha and beta based on theta
        [demand_par,demand_var] = inventory_pars2demand_pars(time_t0_train,xi_vector_train,demand_train,pars_next);
        % weight
        weight=[1/demand_var;1/inventory_var];
        % record parameter estimates for the current iteration
        history(:,iter) = pars_next;
        pars_prev=pars_next;
    end
    pars=pars_next;
%     disp(pars_true)
%     disp(pars0)
%     disp(pars)
    lambda_estimate_ini = pars(1);
    eta_estimate_ini = pars(2);
    rho_estimate_ini = pars(3);
    demand_estimate_ini = demand_par0;
    pars_ini(1:4,j) = [pars0;demand_estimate_ini];
    %
    lambda_estimate = pars(1);
    eta_estimate = pars(2);
    rho_estimate = pars(3);
    demand_estimate = demand_par;
    pars_est(1:4,j) = [pars;demand_estimate];
    save(".\data\parameter.mat","pars_true","pars_ini","pars_est","time0")
    %% fitted level
    time_fit_ini = {};
    time_t0_fit_ini{i} = {};
    demand_fit_ini = {};
    level_diff_fit_ini = {};
    level_fit_ini = {};
    level_fit_Q_ini = {};
    %
    time_fit_est = {};
    time_t0_fit_est{i} = {};
    demand_fit_est = {};
    level_diff_fit_est = {};
    level_fit_est = {};
    level_fit_est_Q = {};
    %
    for i = 1:m
        [time_fit_ini_i,demand_fit_ini_i,level_diff_fit_ini_i,level_fit_ini_i] = inventory_level(demand_estimate_ini,lambda_estimate_ini,eta_estimate_ini,rho_estimate_ini,xi_vector(i),time0,delta_t,Q_vector(i));
        time_fit_ini{i}=time_fit_ini_i;
        time_t0_fit_ini{i} = [time0;time_fit_ini_i];
        demand_fit_ini{i} = demand_fit_ini_i;
        level_diff_fit_ini{i}=level_diff_fit_ini_i;
        level_fit_ini{i}=level_fit_ini_i;
        level_fit_ini_Q{i}=[Q_vector(i);level_fit_ini_i];
        %
        common_value_ini = intersect(time_fit_ini{i},time_true{i});
        indice_time_fit_ini = ismember(time_fit_ini{i}, common_value_ini);
        indice_time_true_ini = ismember(time_true{i}, common_value_ini);
        demand_error_ini_i = abs(demand_fit_ini{i}(indice_time_fit_ini) - demand_true{i}(indice_time_true_ini));
        demand_error_ini(i,j) = mean(demand_error_ini_i);
        level_diff_error_ini_i = abs(level_diff_fit_ini{i}(indice_time_fit_ini) - level_diff_true{i}(indice_time_true_ini));
        level_diff_error_ini(i,j) = mean(level_diff_error_ini_i);
        level_error_ini_i = abs(level_fit_ini{i}(indice_time_fit_ini) - level_true{i}(indice_time_true_ini));
        level_error_ini(i,j) = mean(level_error_ini_i);
        %%
        [time_fit_i,demand_fit_i,level_diff_fit_i,level_fit_i] = inventory_level(demand_estimate,lambda_estimate,eta_estimate,rho_estimate,xi_vector(i),time0,delta_t,Q_vector(i));
        time_fit_est{i}=time_fit_i;
        time_t0_fit_est{i} = [time0;time_fit_i];
        demand_fit_est{i} = demand_fit_i;
        level_diff_fit_est{i}=level_diff_fit_i;
        level_fit_est{i}=level_fit_i;
        level_fit_est_Q{i}=[Q_vector(i);level_fit_i];
        %
        common_value_est = intersect(time_fit_est{i},time_true{i});
        indice_time_fit_est = ismember(time_fit_est{i}, common_value_est);
        indice_time_true_est = ismember(time_true{i}, common_value_est);
        demand_error_est_i = abs(demand_fit_est{i}(indice_time_fit_est) - demand_true{i}(indice_time_true_est));
        demand_error_est (i,j) = mean(demand_error_est_i);
        level_diff_error_est_i = abs(level_diff_fit_est{i}(indice_time_fit_est) - level_diff_true{i}(indice_time_true_est));
        level_diff_error_est(i,j) = mean(level_diff_error_est_i);
        level_error_i = abs(level_fit_est{i}(indice_time_fit_est) - level_true{i}(indice_time_true_est));
        level_error_est(i,j) = mean(level_error_i);
        %
    end
end
% boxplot
xtitle = ["(a) lambda","(b) eta","(c) rho","(d) demand"];
labels = ["initial ests", "Final ests"];
fig_pars = figure('unit','centimeters','position',[10,10,40,10],'PaperPosition',[0, 0, 40,10],'PaperSize',[40,10]);
for i = 1:4
    ax = subplot(1,4,i);
    par_compare = [pars_ini(i,:)', pars_est(i,:)'];
    bp = boxplot(par_compare,'Labels',labels);
    yline(pars_true(i),'r--','LineWidth',1.5);
    set(findobj(gca,'Type','Line'),'LineWidth',1.5);
    set(gca,'FontName','Book Antiqua','FontSize',12)
    title(ax, xtitle(i), 'FontSize',14);
end
fig_error = figure('unit','centimeters', ...
    'position',[10,5,40,20], ...
    'PaperPosition',[0, 0, 40,20], ...
    'PaperSize',[40,20]);
%% ===================== 1. demand train =====================
subplot(2,3,1);
demand_error_ini_train = demand_error_ini(1:6,:);
demand_error_est_train = demand_error_est(1:6,:);
demand_error_ini_train = demand_error_ini_train(:);
demand_error_est_train = demand_error_est_train(:);
demand_error_train = [demand_error_ini_train,demand_error_est_train];
boxplot(demand_error_train,'Labels',labels);
ylim([0,10]);

set(findobj(gca,'Type','Line'),'LineWidth',1.5);
set(gca,'FontName','Book Antiqua','FontSize',12);
title("(a) demand train error",'FontSize',14);
%% ===================== 2. level diff train =====================
subplot(2,3,2);
level_diff_error_ini_train = level_diff_error_ini(1:6,:);
level_diff_error_est_train = level_diff_error_est(1:6,:);
level_diff_error_ini_train = level_diff_error_ini_train(:);
level_diff_error_est_train = level_diff_error_est_train(:);
level_diff_error_train = [level_diff_error_ini_train,level_diff_error_est_train];
boxplot(level_diff_error_train,'Labels',labels);
ylim([0,10]);
set(findobj(gca,'Type','Line'),'LineWidth',1.5);
set(gca,'FontName','Book Antiqua','FontSize',12);
title("(b) level change train error",'FontSize',14);
%% ===================== 3. level train =====================
subplot(2,3,3);
level_error_ini_train = level_error_ini(1:6,:);
level_error_est_train = level_error_est(1:6,:);
level_error_ini_train = level_error_ini_train(:);
level_error_est_train = level_error_est_train(:);
level_error_train = [level_error_ini_train,level_error_est_train];
boxplot(level_error_train,'Labels',labels);
ylim([0,20]);
set(findobj(gca,'Type','Line'),'LineWidth',1.5);
set(gca,'FontName','Book Antiqua','FontSize',12);
title("(c) level train error",'FontSize',14);
%% ===================== 4. demand test =====================
subplot(2,3,4);
demand_error_ini_test = demand_error_ini(7:end,:);
demand_error_est_test = demand_error_est(7:end,:);
demand_error_ini_test = demand_error_ini_test(:);
demand_error_est_test = demand_error_est_test(:);
demand_error_test = [demand_error_ini_test,demand_error_est_test];
boxplot(demand_error_test,'Labels',labels);
ylim([0,10]);
set(findobj(gca,'Type','Line'),'LineWidth',1.5);
set(gca,'FontName','Book Antiqua','FontSize',12);
title("(d) demand test error",'FontSize',14);
%% ===================== 5. level diff test =====================
subplot(2,3,5);
level_diff_error_ini_test = level_diff_error_ini(7:end,:);
level_diff_error_est_test = level_diff_error_est(7:end,:);
level_diff_error_ini_test = level_diff_error_ini_test(:);
level_diff_error_est_test = level_diff_error_est_test(:);
level_diff_error_test = [level_diff_error_ini_test,level_diff_error_est_test];
boxplot(level_diff_error_test,'Labels',labels);
ylim([0,10]);
set(findobj(gca,'Type','Line'),'LineWidth',1.5);
set(gca,'FontName','Book Antiqua','FontSize',12);
title("(e) level change test error",'FontSize',14);
%% ===================== 6. level test =====================
subplot(2,3,6);
level_error_ini_test = level_error_ini(7:end,:);
level_error_est_test = level_error_est(7:end,:);
level_error_ini_test = level_error_ini_test(:);
level_error_est_test = level_error_est_test(:);
level_error_test = [level_error_ini_test,level_error_est_test];
boxplot(level_error_test,'Labels',labels);
ylim([0,20]);
set(findobj(gca,'Type','Line'),'LineWidth',1.5);
set(gca,'FontName','Book Antiqua','FontSize',12);
title("(f) level test error",'FontSize',14);
%
savefig(fig_pars,'.\figure\simulation_pars.fig')
exportgraphics(fig_pars,'.\figure\simulation_pars.pdf')
savefig(fig_error,'.\figure\simulation_error.fig')
exportgraphics(fig_error,'.\figure\simulation_error.pdf')



