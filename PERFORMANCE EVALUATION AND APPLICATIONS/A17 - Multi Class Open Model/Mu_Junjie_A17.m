%MU JUNJIE ASSIGNMENT17
clear all;
L = [2; 3; 2.5] ;
D = [10, 12;
      4, 3;
      6, 6] / 60;
%Utilization of the two stations
Uck =  D.*[L L];
U = sum(Uck);
%Residence time
Rck = D ./ (1 - ones(size(D,1),1)*U);
%Average number of jobs in the system for each type of product
Nck = L .* Rck;
N = sum(Nck,2);
%Average system response time per product type
Rk = sum(Rck,2)*60;
%class-independent average system response time
R = sum(L/60 .* (sum(Rck,2)*60)) / sum(L/60);
fprintf(1, "Utilization of the two stations: %f, %f\n",  U(1),U(2));
fprintf(1, "Average number of jobs in the system for each type of product: %f, %f, %f\n",  N(1),N(2),N(3));
fprintf(1, "Average system response time per product type: %f, %f, %f\n",  Rk(1),Rk(2),Rk(3));
fprintf(1, "Class-independent Aaverage system response time: %f\n",  R);