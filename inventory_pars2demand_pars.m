function [demand_pars,res_var] = inventory_pars2demand_pars(time0,time_train,xi_vector_train,p_vector_train,demand_train,inventory_pars)
% calculate the initial values of alpha and beta based on lambda
% input parameter:
% time0: the time of order arrival
% time_train: the sample time
% xi_vector_train: preservation investment
% p_vector_train: the price vector
% demand_train: the demand
% inventory pars: lambda, rho
% output:
% demand pars: alpha, beta
% res_var: the residual variance of demand regression equation


lambda0=inventory_pars(1);
rho=inventory_pars(2);
% number of orders
cell_length=length(time_train);
E_vector=[];
P_vector=[];
demand=[];
for i = 1:cell_length
    time_i=[time0;time_train{i}];
    demand_i=demand_train{i};
    time_i_j=time_i(2:end);
    time_i_j1=time_i(1:end-1);
    lambda = lambda0*exp(-rho*xi_vector_train(i));
    E_i=lambda\(exp(-lambda*(time_i_j1-time0))-exp(-lambda*(time_i_j-time0)));
    P_i=-p_vector_train(i)*E_i;
    % prepare estimate matrix
    E_vector=[E_vector;E_i];
    P_vector=[P_vector;P_i];
    demand=[demand;demand_i];
end
% least squares
H_lambda=[E_vector,P_vector];
demand_pars=(H_lambda'*H_lambda)\H_lambda'*demand;
% residual
residual =  demand-H_lambda*demand_pars;
res_var = length(residual)\sum(residual.^2);

end

