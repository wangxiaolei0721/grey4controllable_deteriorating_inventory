function res_var = inventory_pars2var(xi_vector_train,demand_train,Q_vector_train,level_train,inventory_pars)
% calculate the residual variance of inventory regression equation based on theta
% input parameter:
% xi_vector_train: preservation investment
% demand_train: the demand
% Q_vector_train: order quantity vector
% level_train: the inventory levels
% inventory pars: lambda, rho
% output parameter:
% res_var: residual variance


% number of orders
cell_length=length(demand_train);
% deterioration pars
lambda0=inventory_pars(1);
rho=inventory_pars(2);
%
Y_inventory=[];
Y_inventory_fit=[];
for i = 1:cell_length
    Q=Q_vector_train(i);
    level=level_train{i};
    Level_Q=[Q;level];
    I1_cumsum=cumsum(Level_Q(1:end-1,1));
    I2_cumsum=cumsum(Level_Q(2:end,1));
    Level_cum=1/2*I1_cumsum+1/2*I2_cumsum;  % accumulative series
    % inventory equation
    demand_i=demand_train{i};
    g_integral=cumsum(demand_i);
    y_inventory=level-Q+g_integral;
    lambda = lambda0*exp(-rho*xi_vector_train(i));
    y_inventory_fit=-lambda*Level_cum;
    Y_inventory=[Y_inventory;y_inventory];
    % 
    Y_inventory_fit=[Y_inventory_fit;y_inventory_fit];
end
% inventory equation error
residual = Y_inventory-Y_inventory_fit;
res_var = length(residual)\sum(residual.^2);


