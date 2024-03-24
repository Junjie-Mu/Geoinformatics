%MU JUNJIE ASSIGNMENT10
clear all

lambda = 40;
D = 16 /1000 ;
t = 0.5;
k = 90;

rho = lambda * D;
U = rho ; 
%Probability of having exactly one job
P_1 = (1-rho) * rho;
%Probability of having less than 10 jobs 
P_less10 = 1 - rho^(10+1);
%Average queue length
aN = rho / (1-rho);
%Average response time
aR = D / (1-rho);
%Probability the response time > 0.5
P_Rgt = exp (-t/aR);
%90 precentile of the response time distribution
Theta_k = - log(1-k/100)*aR;

fprintf(1, "M/M/1\n");
fprintf(1, " The utilization of the system: %g\n",  U);
fprintf(1, " The probability of having exactly one job in the system: %g\n", P_1 );
fprintf(1, " The probability of having less than 10 jobs in the system: %g\n",  P_less10);
fprintf(1, " The average queue length (job not in service): %g\n",  aN);
fprintf(1, " The average response time: %g\n",  aR);
fprintf(1, " The probability that the response time is greater than 0.5 s: %g\n",  P_Rgt);
fprintf(1, " The 90 percentile of the response time distribution: %g\n\n",  Theta_k);



lambda2 = 90;
D2 = 16 / 1000;
c = 2;
rho2 = lambda2 * D2 / c;

U2 = lambda2 * D2;
aU2 = rho2;
P2_1 = 2 * (1-rho2) * rho2 / (1+rho2);
P2_less10 = 1 - 2*rho2^(10+1)/(1+rho2);
aN2 = 2*rho2/(1-rho2^2);
aR2 = D / (1-rho2^2);


fprintf(1, "M/M/2\n");
fprintf(1, " The total and average utilizations of the system: %g，%g\n",U2, aU2);
fprintf(1, " The probability of having exactly one job in the system: %g\n", P2_1 );
fprintf(1, " The probability of having less than 10 jobs in the system: %g\n",  P2_less10);
fprintf(1, " The average queue length (job not in service): %g\n",  aN2);
fprintf(1, " The average response time: %g\n",  aR2);