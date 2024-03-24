%MU JUNJIE ASSIGNMENT15
clear all;
N = 80;
S = [40000,50,2,80,80,120]*10^(-3);
l = [1, 0, 0, 0, 0 ,0];
p0 = [  0, 1,   0,   0,   0,   0;
        0, 0, 0.4, 0.5,   0,   0;
        0, 0,   0,   0, 0.6, 0.4;
        0, 1,   0,   0,   0,   0;
        0, 1,   0,   0,   0,   0;
        0, 1,   0,   0,   0,   0];
V = l * inv(eye(6) - p0);
D = V' .* S';
Z = 40;

Nk = zeros(6,N);
R = D;

for i = 1: N-1
    X = i/(sum(R(:,i))+Z);
    Nk(:,i) = X*R(:,i);
    R(:,i+1) = (1+Nk(:,i)) .* D;
end
%The throughput of the system (X)
X = N/(sum(R(:,N))+Z);
%The average system response time (R)
sR = sum(R(:,N));
%The utilization of the [3] Application Server, [4] DBMS, [5] Disk 1, and [6] Disk 2.
U = X * D;
fprintf(1, "System Throughput: %f\n",  X);
fprintf(1, "Average system response time: %f\n",  sR); 
fprintf(1, "Utilization of Application Server: %f, " + ...
    "DBMS: %f, Disk 1:%f, Disk 2:%f\n",  U(3,1), U(4,1), U(5,1), U(6,1));


