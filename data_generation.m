% clear data and figure
clc;
clear;
close all;
%% model setting
% demand parameters
alpha=60;
beta=6;
% deterioration parameters
lambda0=0.1;
rho=5;
% order cycles
m=10;
% random seed
rng(10)
% the order quantity
Q_vector=300+randi([0,60],m,1);
% the sales price
% the random scalar obtained from the uniform distribution of interval (0,1)
p_vector=4+2*rand(m,1);
d=alpha-beta*p_vector;
% the preservation cost
xi_vector = 0.2*rand(m,1); % 0.5*rand(m,1)
lambda=lambda0*exp(-rho*xi_vector);
%% simulation parameters
% the time of order arrival
time0=0;
% the time resolution
delta_t=1;
% standard deviation of error in demand regression equation
std_dev=2;
%% initialization of data storage
time_true = {};
time_true_t0 = {};
demand_true = {};
level_true = {};
time_simu = {};
time_simu_t0 = {};
demand_simu = {};
level_simu = {};
% generate the inventory levels
for i = 1:m
    % true level
    [time_true_i,demand_true_i,level_true_i] = inventory_level(alpha,beta,p_vector(i),lambda0,rho,xi_vector(i),time0,delta_t,Q_vector(i));
    time_true{i}=time_true_i;
    time_true_t0{i} = [time0;time_true_i];
    demand_true{i} = demand_true_i; % [alpha-beta*p_vector(i);demand_true_i]
    level_true{i}=level_true_i;
    % simulated level
    [time_simu_i,demand_simu_i,level_simu_i] = inventory_level_simulation(alpha,beta,p_vector(i),std_dev,lambda0,rho,xi_vector(i),time0,delta_t,Q_vector(i));
    time_simu{i} = time_simu_i;
    demand_simu{i} = demand_simu_i;
    time_simu_t0{i} = [time0;time_simu_i];
    level_simu{i}=level_simu_i;
end
%% plot
% demand plot
fdemand=figure('unit','centimeters','position',[5,5,40,20],'PaperPosition',[5,5,40,20],'PaperSize',[40,20]);
tiledlayout(2,m/2,'Padding','Compact');
% plot time vs level
finvertory=figure('unit','centimeters','position',[5,5,40,20],'PaperPosition',[5,5,40,20],'PaperSize',[40,20]);
tiledlayout(2,m/2,'Padding','Compact');
for i = 1:m
    figure(fdemand)
    nexttile
    plot(time_true{i},demand_true{i},'LineWidth',1)
    hold on
    plot(time_simu{i},demand_simu{i},'LineWidth',1)
    xlabel({'Day'},'FontSize',12)
    xlim([min(time_simu{i}),max(time_simu{i})])
    ylabel(['Demand'],'FontSize',12)
    title(strcat("(",char(96 + i),") The ", num2str(i),"th ordering cycle"),'FontSize',14)
    set(gca,'FontName','Book Antiqua','FontSize',10)
    if i==10
        legend(["Actual demand","Simulated demand"],'location','northeast','FontSize',10,'NumColumns',1)
    end
    figure(finvertory)
    nexttile
    plot(time_true{i},level_true{i},'LineWidth',1)
    hold on
    plot(time_simu{i},level_simu{i},'LineWidth',1)
    xlabel({'Day'},'FontSize',12)
    xlim([min(time_simu{i}),max(time_simu{i})])
    ylabel(['Inventory level'],'FontSize',12)
    title(strcat("(",char(96 + i),") The ", num2str(i),"th ordering cycle"),'FontSize',14)
    set(gca,'FontName','Book Antiqua','FontSize',10)
    if i==10
        legend(["Actual inventory level","Simulated inventory level"],'location','northeast','FontSize',8,'NumColumns',1) % ,"Fitted inventory level"
    end
end
%% save data
save(".\data\external_settings.mat","time0","delta_t","m","p_vector","xi_vector","Q_vector")
save(".\data\simulated_data.mat","time_simu","demand_simu","level_simu")
save(".\data\true_data.mat","time_true","time_true_t0","demand_true","level_true")

