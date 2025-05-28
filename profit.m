function profit = profit(alpha,beta,p,lambda0,rho,xi,c,h,K,T)
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
par1=(alpha-beta*p).*(p.*lambda+h+xi)./(lambda.^2);
par2=(alpha-beta*p).*(c*lambda+h+xi)./lambda;
profit=par1.*(1-exp(-lambda.*T))./T-par2-K./T;

end
