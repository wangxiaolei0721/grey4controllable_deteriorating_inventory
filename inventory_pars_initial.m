function inventory_pars = inventory_pars_initial(time_t0_train,xi_train,demand_train,level_diff_train,level_train)
% calculate the initial value of lambda according to the inventory regression equation
% input parameter:
% time_t0_train: the sample time, including time0
% xi_train: preservation investment
% demand_train: the demand
% level_diff_train: simulated inventory changes
% level_train: the inventory levels
% output parameter:
% inventory pars: lambda, eta, rho


% number of orders
cell_length=length(time_t0_train);
time_1f1 = [];
delta_time=[];
xi=[];
inventory_diff=[];
demand=[];
inventory=[];
for i = 1:cell_length
    time_t0_i= time_t0_train{i};  
    time_1f1_i =0.5*time_t0_i(1:(end-1))+0.5*time_t0_i(2:end);
    time_i_diff=diff(time_t0_i);
    xi_i=repmat(xi_train(i),length(time_i_diff),1);
    demand_i=demand_train{i};
    level_diff_i=level_diff_train{i};
    level_i = level_train{i};
    % inventory equation
    time_1f1 = [time_1f1;time_1f1_i];
    delta_time=[delta_time;time_i_diff];
    xi=[xi;xi_i];
    inventory_diff=[inventory_diff;level_diff_i];
    demand=[demand;demand_i];
    inventory_i=0.5*(level_i(2:end)+level_i(1:(end-1)));
    inventory=[inventory;inventory_i];
end
% the initial value of lambda
Y_inventory = log(-inventory_diff-demand) - log(inventory.*delta_time);
H_inventory = [ones(length(delta_time),1),log(time_1f1),-xi];
eq_pars = (H_inventory'*H_inventory)\H_inventory'*Y_inventory;
inventory_pars=[exp(eq_pars(1))/(eq_pars(2)+1);eq_pars(2)+1;eq_pars(3)];

end




