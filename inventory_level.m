function [time,demand,level] = inventory_level(alpha,beta,p,lambda0,rho,xi,time0,delta_t,Q)
% generate inventory levels and inventory changes
% input parameter:
% alpha: basic demand
% beta: price sensitivity coefficient
% p: price
% lambda0: natural deteriorating rate
% rho: parameters in the reduced deterioration rate ratio function
% xi: preservation cost
% time0: the time of order arrival
% delta_t: the time resolution
% Q: the order quantity
% output parameter
% time: sampling time
% demand: demand quantity
% level: inventory level



% calculate order cycle based on order quantity
T=Q/(alpha-beta*p);
% the moment when the inventory drops to 0
tT=time0+T;
% generate sampling time with time resolution as the step size
time=[(time0+delta_t):delta_t:tT]';
% demand quantity
lambda=lambda0*exp(-rho*xi);
demand=(alpha-beta*p)*exp(-lambda*(time-time0));
% inventory levels
level=(alpha-beta*p)*exp(-lambda*(time-time0)).*(tT-time);


end

