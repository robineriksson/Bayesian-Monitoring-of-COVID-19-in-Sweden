function [F_l,H_l,K_l,betaIdx_l] = lumpIAEC19(F,H,K,rates,n)
%LUMPIAEC19 Lump I,A,E in Covid-19 model.
%   [F_l,H_l,K_l,betaIdx_l] = LUMPIAEC19(F,H,K,betaIdx,rates)
%   produces system matrices where the states I,A,E have been lumped to
%   one, such that steady state rates and viral shedding out of the lumped
%   compartment coincide with those of the full model. 

% H. Runvik 2021-01-10

if nargin<5
    n=1;
end
alpha=F(2,3);
beta=F(1,3);
gamma=rates.gammaA(n)*(1-F(1,2)/rates.gammaA(n));
delta=F(1,2);
theta=rates.gammaI(n);

theta_l=theta*(beta*(gamma+delta)+delta*alpha)/(theta*alpha+theta*(gamma+delta)+beta*(gamma+delta)+alpha*delta);
gamma_l=theta*gamma*alpha/(theta*alpha+theta*(gamma+delta)+beta*(gamma+delta)+alpha*delta);
A_rel=theta*alpha/(theta*alpha+theta*(gamma+delta)+beta*(gamma+delta)+alpha*delta);
E_rel=theta*(gamma+delta)/(theta*alpha+theta*(gamma+delta)+beta*(gamma+delta)+alpha*delta);
I_rel=(beta*(gamma+delta)+delta*alpha)/(theta*alpha+theta*(gamma+delta)+beta*(gamma+delta)+alpha*delta);

F_l=zeros(5,5);
F_l(2:end,2:end) = F(4:end,4:end);
F_l(1,2) = F(3,4);
F_l(1,1) = 1-theta_l-gamma_l;
F_l(2,1) = I_rel*F(4,1)+A_rel*F(4,2)+E_rel*F(4,3);
F_l(3,1) = I_rel*F(5,1);
F_l(5,1) = I_rel*F(7,1);

H_l = H(:,3:end);
K_l = K(3:end,:);
betaIdx_l=false(5,5);
betaIdx_l(1,2)=true;


