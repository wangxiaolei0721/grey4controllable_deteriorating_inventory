function [demand_par,res_var] = inventory_pars2demand_pars(time_t0_train,xi_train,demand_train,inventory_pars)
% calculate the initial values of alpha and beta based on lambda
% input parameter:
% time_t0_train: the sample time, including time0
% xi_train: preservation investment
% demand_train: the demand
% inventory pars: lambda, eta, rho
% output:
% demand pars: d
% res_var: the residual variance of demand regression equation


% pars
lambda = inventory_pars(1);
eta = inventory_pars(2);
rho = inventory_pars(3);
% number of orders
cell_length=length(time_t0_train);
E_vector=[];
P_vector=[];
demand=[];
for i = 1:cell_length
    time_t0_i=time_t0_train{i};
    demand_i=demand_train{i};
    time_1f1_i =0.5*(time_t0_i(1:(end-1))+time_t0_i(2:end));
    E_i=exp(-lambda*exp(-rho*xi_train(i)).*(time_1f1_i).^eta);
    % prepare estimate matrix
    E_vector=[E_vector;E_i];
    demand=[demand;demand_i];
end
% least squares
H_par=E_vector;
demand_par=(H_par'*H_par)\H_par'*demand;
% residual
residual =  demand-H_par*demand_par;
res_var = length(residual)\sum(residual.^2);

end

