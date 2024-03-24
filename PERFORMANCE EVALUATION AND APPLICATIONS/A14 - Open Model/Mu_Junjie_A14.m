%MU JUNJIE ASSIGNMENT114
clear all;
%Scenario 1
lambdaIn = [3, 2, 0, 0];
lambda0 = sum(lambdaIn);
P = [   0,  0.8,    0,     0;
        0,    0,  0.3,   0.5;
        0,    1,    0,     0;
        0,    1,    0,    0];
l = lambdaIn/lambda0;
%Visits to the four stations
V = l * inv(eye(4) - P);
S = [2000,30,100,80]*10^(-3);
%Demand of the four stations
D = V .* S;
U = lambda0 * D;
%Throughput of the system
X = lambda0;
for i = 1: 4
    Ni(1,i) = U(i) / (1-U(i));
    Ri(1,i) = D(i) / (1-U(i));
end
%Average number of jobs in the system
N = sum(Ni);
%Average system response time
R = sum(Ri);

fprintf(1, "Throughput of the system:%f\n", X);
fprintf(1, "Average number of jobs in the system:%f\n", N);
fprintf(1, "Average system response time:%f\n", R);