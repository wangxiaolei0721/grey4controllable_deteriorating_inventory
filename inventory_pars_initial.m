function inventory_pars = inventory_pars_initial(xi_vector_train,demand_train,Q_vector_train,level_train)
% calculate the initial value of lambda according to the inventory regression equation
% input parameter:
% xi_vector_train: preservation investment
% demand_train: the demand
% Q_vector_train: order quantity vector
% level_train: the inventory levels
% output parameter：
% res_var: residual variance


cell_length=length(demand_train);
xi=[];
Y_inventory=[];
for i = 1:cell_length
    Q=Q_vector_train(i);
    Level=level_train{i};
    Level_Q=[Q;Level];
    I1_cumsum=cumsum(Level_Q(1:end-1,1));
    I2_cumsum=cumsum(Level_Q(2:end,1));
    Level_cum=1/2*I1_cumsum+1/2*I2_cumsum;  % accumulative series
    % estimate parameter
    demand=demand_train{i};
    g_integral=cumsum(demand);
    y_inventory=log(-(Level-Q+g_integral))-log(Level_cum);
    Y_inventory=[Y_inventory;y_inventory];
    % 
    xi_i=-repmat(xi_vector_train(i),length(Level),1);
    xi=[xi;xi_i];
end

H_inventory=[ones(length(xi),1),xi];
eq_par = (H_inventory'*H_inventory)\H_inventory'*Y_inventory;
inventory_pars=[exp(eq_par(1));eq_par(2)];
%
% residual =  Y_inventory-H_inventory*eq_par;
% inventory_res_var = length(residual)\sum(residual.^2);


end




