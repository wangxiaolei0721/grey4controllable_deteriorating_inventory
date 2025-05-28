function profit = profit_appro(alpha,beta,p,lambda0,rho,xi,c,h,K,T)
% profit function
% input parameter:
% alpha: basic demand
% beta: price sensitivity coefficient
% p: price
% lambda0: natural deteriorating rate
% rho: parameters in the reduced deterioration rate ratio function
% xi: preservation cost
% c: production cost
% h: holding cost per unit per unit of time
% K: ordering cost per cycle
% T: order cycles
% output parameter:
% profit: total profits at T


lambda=lambda0*exp(-rho*xi);
profit=(alpha-beta*p).*(p-c)-(alpha-beta*p).*(p.*lambda+h+xi).*T/2-K./T;


end
