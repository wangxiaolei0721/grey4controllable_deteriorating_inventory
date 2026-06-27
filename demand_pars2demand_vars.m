function res_var = demand_pars2demand_vars(time_t0_test,xi_test,demand_test,inventory_pars,demand_pars)
% calculate the initial values of alpha and beta based on lambda
% input parameter:
% time_t0_test: the sample time, including time0
% xi_vector_test: preservation investment
% demand_test: the demand data
% inventory pars: lambda, eta, rho
% demand par: d
% output:
% demand_var: the residual variance of demand regression equation

% pars
lambda = inventory_pars(1);
eta = inventory_pars(2);
rho = inventory_pars(3);
% number of orders
cell_length=length(time_t0_test);
E_vector=[];
P_vector=[];
demand=[];
for i = 1:cell_length
    time_t0_i=time_t0_test{i};
    demand_i=demand_test{i};
    time_1f1_i =0.5*time_t0_i(1:(end-1))+0.5*time_t0_i(2:end);
    E_i=exp(-lambda*exp(-rho*xi_test(i)).*(time_1f1_i).^eta);
    % prepare estimate matrix
    E_vector=[E_vector;E_i];
    demand=[demand;demand_i];
end
% least squares
H_par=E_vector;
% residual
residual =  demand-H_par*demand_pars;
res_var = length(residual)\sum(residual.^2);

end

