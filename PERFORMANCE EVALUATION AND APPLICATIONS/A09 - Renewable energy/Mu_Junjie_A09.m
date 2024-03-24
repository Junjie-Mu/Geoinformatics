%MU JUNJIE ASSIGNMENT9
clear all

MTTF1 = 18;
MTTF2 = 3;
MTTF3 = 8;
MTTR1 = 2;
MTTR2 = 2;
MTTR3 = 2;

l1 = 1/MTTF1;
l2 = 1/MTTF2;
l3 = 1/MTTF3;
m1 = 1/MTTR1;
m2 = 1/MTTR2;
m3 = 1/MTTR3;

Q = [ -l1-l2-l3, l1 , l2  , l3;
      m1       , -m1, 0   , 0;
      m2       , 0  , -m2 , 0;
      m3       , 0  , 0   , -m3]

p0 = [0.25, 0.25 , 0.25 ,0.25];
alpha1 = [0.1, 12, 12, 12]
alpha2 = [0, 1/2, 1/4, 1/4]
Tr1 = [ 0,0,0,0;
        1/2,0,0,0;
        1/4,0,0,0;
        1/4,0,0,0]
[t, Sol] = ode45(@(t,x) Q'*x, [0 1440], p0');

plot(t, Sol, "-");
legend("SCAN", "NIGHTSLEEP", "SUNNYSLEEP", "CLOUDYSLEEP");


power = Sol * alpha1';
utilization = Sol * alpha2';
throughput = Sol * Tr1';


fprintf("Average Power Consumption: %f\n", mean(power));
fprintf("Average Utilization: %f\n", sum(utilization)/1440);
fprintf("Throughput: %f", mean(sum(throughput)));

Sol(end,:)
[Sol(end,:) * alpha1', max(Sol * alpha1')]
[Sol(end,:) * alpha2', max(Sol * alpha2')]
[Sol(end,:) * Tr1', min(Sol * Tr1')]
