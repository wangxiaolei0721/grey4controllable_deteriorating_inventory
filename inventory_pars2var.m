function res_var = inventory_pars2var(time_t0_train,xi_train,demand_train,level_diff_train,level_train,inventory_pars)
% calculate the residual variance of inventory regression equation based on theta
% input parameter:
% time_t0_train: the sample time, including time0
% xi_train: preservation investment
% demand_train: the demand
% level_diff_train: simulated inventory changes
% level_train: the inventory levels
% inventory pars: lambda, eta, rho
% output parameter:
% res_var: residual variance


% pars
lambda = inventory_pars(1);
eta = inventory_pars(2);
rho = inventory_pars(3);
% number of orders
cell_length=length(time_t0_train);
time_1f1 = [];
delta_time=[];
inventory_diff=[];
demand=[];
det_qua=[];
for i = 1:cell_length
    time_t0_i= time_t0_train{i};  
    time_1f1_i =0.5*time_t0_i(1:(end-1))+0.5*time_t0_i(2:end);
    time_i_diff=diff(time_t0_i);
    demand_i=demand_train{i};
    level_diff_i=level_diff_train{i};
    level_i = level_train{i};
    % inventory equation
    delta_time=[delta_time;time_i_diff];
    inventory_diff=[inventory_diff;level_diff_i];
    demand=[demand;demand_i];
    inventory_i=0.5*(level_i(2:end)+level_i(1:(end-1))); % 0.5*(level_i(2:end)+level_i(1:end-1));
    theta_t_i = lambda*eta*exp(-rho*xi_train(i))*(time_1f1_i).^(eta-1);
    det_qua_i=-theta_t_i.*inventory_i;
    det_qua=[det_qua;det_qua_i];
end
% inventory equation error using nonlinear equation
Y_inventory = inventory_diff+demand;
% inventory residual
residual =  Y_inventory-det_qua.*delta_time;
res_var = (length(residual)-1)\sum(residual.^2);
