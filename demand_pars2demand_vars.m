function res_var = demand_pars2demand_vars(time0,time_test,xi_vector_test,p_vector_test,demand_test,inventory_pars,demand_pars)
% calculate the initial values of alpha and beta based on lambda
% input parameter:
% time0: the time of order arrival
% time_train: the sample time
% xi_vector_test: preservation investment
% demand_train: the demand data
% p_vector_train: the price vector
% inventory pars: lambda, rho
% demand pars: alpha, beta
% output:
% demand_var: the residual variance of demand regression equation

lambda0=inventory_pars(1);
rho=inventory_pars(2);
% number of orders
cell_length=length(time_test);
E_vector=[];
P_vector=[];
demand=[];
for i = 1:cell_length
    time_i=[time0;time_test{i}];
    demand_i=demand_test{i};
    time_i_j=time_i(2:end);
    time_i_j1=time_i(1:end-1);
    lambda = lambda0*exp(-rho*xi_vector_test(i));
    E_i=lambda\(exp(-lambda*(time_i_j1-time0))-exp(-lambda*(time_i_j-time0)));
    P_i=-p_vector_test(i)*E_i;
    % prepare estimate matrix
    E_vector=[E_vector;E_i];
    P_vector=[P_vector;P_i];
    demand=[demand;demand_i];
end
% least squares
H_lambda=[E_vector,P_vector];
% residual
residual =  demand-H_lambda*demand_pars;
res_var = length(residual)\sum(residual.^2);

end

