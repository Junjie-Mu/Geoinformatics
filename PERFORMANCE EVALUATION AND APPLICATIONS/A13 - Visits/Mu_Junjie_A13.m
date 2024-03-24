%MU JUNJIE ASSIGNMENT13
clear all;
%Scenario 1
lambdaIn = [1, 0, 0, 0];
lambda0 = sum(lambdaIn);
P = [   0,    1,    0,     0;
      0.1,    0,  0.3,   0.6;
        0, 0.85,    0,  0.15;
        0, 0.75, 0.25,    0];
Pn = [zeros(4,1), P(:,2:end)];
l = lambdaIn/lambda0;
V = l * inv(eye(4) - Pn);
%Visits to the four stations
S = [10,20,10,3]*10^(-3);
%Demand of the four stations
D = V .* S;
fprintf(1, "Scenario 1\n");
fprintf(1, "[Visits] Terminals:%f, CPU:%f, Disk:%f, RAM:%f\n", V(1,1),V(1,2),V(1,3),V(1,4));
fprintf(1, "[Demand] Terminals:%f, CPU:%f, Disk:%f, RAM:%f\n", D(1,1),D(1,2),D(1,3),D(1,4));

%Scenario 2
lambdaIn_2 = [0.3, 0, 0];
lambda0_2 = sum(lambdaIn_2);
P_2 = [      0,     0.3,   0.6;
         0.8,       0,  0.15;
        0.75,    0.25,    0;];
l_2 = lambdaIn_2/lambda0_2;
V_2 = l_2 * inv(eye(3) - P_2);
%Visits to the three stations
S_2 = [20,10,3]*10^(-3);
%Demand of the three stations
D_2 = V_2 .* S_2;
%Throughput of the three stations
X = 0.3 * V_2;
fprintf(1, "\nScenario 2\n");
fprintf(1, "[Visits] CPU:%f, Disk:%f, RAM:%f\n", V_2(1,1),V_2(1,2),V_2(1,3));
fprintf(1, "[Demand] CPU:%f, Disk:%f, RAM:%f\n", D_2(1,1),D_2(1,2),D_2(1,3));
fprintf(1, "[Throughput] CPU:%f, Disk:%f, RAM:%f\n", X(1,1),X(1,2),X(1,3));