function inventory_pars_res = inventory_pars_res(time0,time_train,xi_vector_train,p_vector_train,Q_vector_train,demand_train,level_train,weight,inventory_pars)
% Iteratively Reweighed Least Squares algorithm
% input parameter:
% time0: the time of order arrival
% time_train: the sample time
% p_vector_train: the price vector
% Q_vector_train: order quantity vector
% demand_train: the demand
% level_train: the inventory levels
% weight_initial: initial weight
% inventory pars: lambda, rho
% output:
% inventory_pars_res: the inventory residual and the demand residual



% number of orders
cell_length=length(time_train);
% deterioration pars
lambda0=inventory_pars(1);
rho=inventory_pars(2);
%
Y_inventory=[];
Y_inventory_fit=[];
% demand equation
E_vector=[];
P_vector=[];
demand=[];
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
    % demand equation
    time_i=[time0;time_train{i}];
    time_i_j=time_i(2:end);
    time_i_j1=time_i(1:end-1);
    E_i=lambda\(exp(-lambda*(time_i_j1-time0))-exp(-lambda*(time_i_j-time0)));
    P_i=-p_vector_train(i)*E_i;
    % prepare estimate matrix
    E_vector=[E_vector;E_i];
    P_vector=[P_vector;P_i];
    demand=[demand;demand_i];
end


% inventory residual
inventory_res =  weight(2)*(Y_inventory-Y_inventory_fit);
% demand equation error
% least squares
H_lambda=[E_vector,P_vector];
demand_par=(H_lambda'*H_lambda)\H_lambda'*demand;
% residual
demand_residual =  weight(1)*(demand-H_lambda*demand_par);
%
inventory_pars_res=[inventory_res;demand_residual];


end

