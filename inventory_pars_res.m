function inventory_pars_res = inventory_pars_res(time_t0_train,xi_train,demand_train,level_diff_train,level_train,weight,inventory_pars)
% Iteratively Reweighed Least Squares algorithm
% input parameter:
% time_t0_train: the sample time, including time0
% xi_train: preservation investment
% demand_train: the demand
% level_diff_train: simulated inventory changes
% level_train: the inventory levels
% weight_initial: initial weight
% inventory pars: lambda, eta, rho
% output:
% inventory_pars_res: the inventory residual and the demand residual

% output parameter:
% inventory pars: lambda, eta, rho


% number of orders
cell_length=length(time_t0_train);
% deterioration pars
lambda = inventory_pars(1);
eta = inventory_pars(2);
rho = inventory_pars(3);
%
inventory_diff =[];
inventory_fit=[];
% demand equation
E_vector=[];
demand=[];
delta_time=[];
for i = 1:cell_length
    time_t0_i= time_t0_train{i};
    time_1f1_i =0.5*(time_t0_i(1:(end-1))+time_t0_i(2:end));
    time_i_diff=diff(time_t0_i);
    level_diff_i=level_diff_train{i};
    level_i = level_train{i};
    % inventory equation
    delta_time=[delta_time;time_i_diff];
    inventory_diff=[inventory_diff;level_diff_i];
    inventory_i=0.5*(level_i(2:end)+level_i(1:(end-1))); % 0.5*(level_i(2:end)+level_i(1:end-1));
    theta_t_i = lambda*eta*exp(-rho*xi_train(i))*(time_1f1_i).^(eta-1);
    inventory_fit_i=-theta_t_i.*inventory_i;
    inventory_fit=[inventory_fit;inventory_fit_i];
    % demand equation
    demand_i=demand_train{i};
    E_i=exp(-lambda*exp(-rho*xi_train(i)).*(time_1f1_i).^eta);
    % prepare estimate matrix
    E_vector=[E_vector;E_i];
    demand=[demand;demand_i];
end


inventory = inventory_diff+demand;
% inventory residual
inventory_res =  weight(2)*(inventory-inventory_fit.*delta_time);
% demand equation error
% least squares
H_par=E_vector;
demand_par=(H_par'*H_par)\H_par'*demand;
% residual
demand_residual =  weight(1)*(demand-H_par*demand_par);
%
inventory_pars_res=[inventory_res;demand_residual];

end

