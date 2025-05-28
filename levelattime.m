function level = levelattime(alpha,beta,p,lambda0,rho,xi,time,time0,Q)
% calculate inventory levels at a specific time
% input parameter:
% alpha: basic demand
% beta: price sensitivity coefficient
% p: price
% lambda0: natural deteriorating rate
% rho: parameters in the reduced deterioration rate ratio function
% xi: preservation cost
% time: sampling time
% time0: the time of order arrival
% delta_t: the time resolution
% Q: the order quantity
% output parameter
% level: inventory level


lambda=lambda0*exp(-rho*xi);
% calculate order cycle based on order quantity
T=Q/(alpha-beta*p);
% the moment when the inventory drops to 0
tT=time0+T;
level=(alpha-beta*p)*exp(-lambda*(time-time0)).*(tT-time);

end
