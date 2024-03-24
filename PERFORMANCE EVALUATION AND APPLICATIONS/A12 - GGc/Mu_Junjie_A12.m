%MU JUNJIE ASSIGNMENT12
clear all;

lambda = 10;
mu1 = 50;
mu2 = 5;
p1 = 0.8;

D = p1/mu1 + (1-p1)/mu2; 
m2 = 2 *(p1/mu1^2 + (1-p1)/mu2^2);
Var = m2 - D^2;
cv2 = Var / D^2;
rho = lambda *D ;
%Utilization of the system
U = lambda * D;
%Average response time
R = D + (rho*D) / (1-rho) * (1+cv2) / 2;
%Average number of jobs
N = rho + rho^2 * (1+cv2) / (2*(1-rho));

fprintf(1, "The utilization of the system: %g\n",  U);
fprintf(1, "The (exact) average response time: %g\n",  R);
fprintf(1, "The (exact) average number of jobs in the system: %g\n",  N);

lamda = 240;
k = 5;
c = 3;
T = k/lamda;
lamda = 1/T;

Da = k / lamda;
ma2 = k / lamda^2;
ca2 = (ma2 -Da^2) / Da^2;
rho = D / (c*T);

%Approxiamte average response time
theta = D/(c*(1-rho))/(1+(1-rho)*(factorial(c)/(c*rho)^c)*(sum((c*rho).^(0:c-1) ./factorial(0:2))));
R = D + (ca2+ cv2)/2 * theta;
%Approxiamte average number of jobs in the system
N = lamda/k *R;

fprintf(1, "\nThe average utilization of the system: %g\n",  rho);
fprintf(1, "The approxiamte average response time: %g\n",  R);
fprintf(1, "The approxiamte average number of jobs in the system: %g\n",  N);
