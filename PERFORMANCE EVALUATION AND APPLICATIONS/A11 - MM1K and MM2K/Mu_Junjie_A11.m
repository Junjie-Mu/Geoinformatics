%MU JUNJIE ASSIGNMENT11
clear all

c = 1;
k = 32;
lambda = 150;
D = 350/1000;

rho = lambda * D;
%Utilization of the system
U = (rho - rho^(k+1)) / ( 1 - rho^(k+1)) ;
%Loss probability
LP = (rho^k - rho^(k+1)) / (1 - rho^(k+1));
%Average number of jobs in the system
N = rho/(1-rho) - (k+1)*rho^(k+1) / (1-rho^(k+1));
%Drop rate
DR = lambda * (rho^k - rho^(k+1))/(1-rho^(k+1));
%Average response time
R = D * (1-(k+1)*rho^k+k*rho^(k+1)) / (1-rho) / (1-rho^k);
%Average time spent in the queue(waiting for service)
aTS = R - D;

fprintf(1, "M/M/1/K\n");
fprintf(1, "Utilization of the system: %f\n",  U);
fprintf(1, "Loss probability: %f\n",  LP);
fprintf(1, "Average number of jobs in the system: %f\n",  N);
fprintf(1, "Drop rate: %f\n",  DR);
fprintf(1, "Average response time: %f\n",  R);
fprintf(1, "Average time spent in the queuee: %f\n",  aTS);

c = 2;
lambda = 250;
rho = lambda * D / c;
P0 = ((c*rho)^c/factorial(c)*(1-rho^(k-c+1)) / (1-rho) + 1 + c*rho) ^(-1);
%Utilization
U2 = P0 / factorial(1) *(c*rho)^1 + 2 * (P0 / factorial(2) *(c*rho)^2);
for i = (c+1) : k
    Pi = P0 * c^c * rho^i / factorial(c);
    U2 = U2 + c*Pi;
end
P0 = ((c*rho)^c/factorial(c)*(1-rho^(k-c+1)) / (1-rho)) ^(-1);
Pk = P0 * c^c * rho^k / factorial(c);
%Loss probability
Pl = Pk;
%Average number of jobs
N = 0;
for i = 1 : k
    if i < c
        Pi = P0 / factorial(i) *(c*rho)^i;
        N = N + i*Pi;
    else
        Pi = P0 * c^c * rho^i / factorial(c);
        N = N + i*Pi;
    end
end
%Drop rate
Dr = lambda * Pk;
%Average response time
R = N / (lambda * (1-Pk));
%Average time spent in the queue
aTS = R - D;

fprintf(1, "\nM/M/2/K\n");
fprintf(1, "Total and average utilization of the channel: %f,%f\n", U2, U2/2);
fprintf(1, "Loss probability: %f\n",  Pl);
fprintf(1, "Average number of jobs in the system: %f\n",  N);
fprintf(1, "Drop rate: %f\n",  Dr);
fprintf(1, "Average response time: %f\n",  R);
fprintf(1, "Average time spent in the queue: %f\n",  aTS);
