clear all;

filename = 'barrier.log';
file = fopen(filename,"r");
AC = cell(0,2);
while ~feof(file)
    line = fgetl(file);
    match = regexp(line, '\[(\d+:\d+:\d+:\d+:\d+:\d+)\]\[(\w+)\]', 'tokens');  
    if ~isempty(match)
        timestamp = match{1}{1}; 
        event = match{1}{2}; 
        parts = sscanf(timestamp, '%d:%d:%d:%d:%d:%d', 6);
        timeInSeconds = parts(2)* 86400 + parts(3) * 3600 + parts(4) * 60 + parts(5) + parts(6)/10;
        if strcmp(event, '_IN_')
            %in_time = timeInSeconds;
            AC{end+1 ,1} =timeInSeconds;
        elseif strcmp(event, '_OUT')
            AC{end+1 ,2} =timeInSeconds;
        end
    end
end
fclose(file);
AC = cell2mat(AC);
%Arrivals
nA = size(AC, 1);
%Completions
nC = size(AC, 1);
T = AC(end, 2) - AC(1, 1);
%Arrival rate
Lambda = nA / T;
%Throughput
X = nC / T;
A = AC(:,1);
C= AC(:,2);
It = A(2: end) - A(1: end -1);
an = C(end,1) - A(end,1);
It(end+1 ,1) = C(end,1) - A(end,1);
%Average inter-arrival time
At = sum(It)  / nA;
Rt = C - A;
W = sum(Rt);

%Average response time
R = W / nC;
%Average number of jobs
N = W / T;

evs = [AC(:,1), ones(nA, 1), zeros(nA, 4); AC(:,2), -ones(nC,1), zeros(nC, 4)];
evs = sortrows(evs, 1);
evs(:,3) = cumsum(evs(:,2));
evs(1:end-1, 4) = evs(2:end,1) - evs(1:end-1,1);
evs(:,5) = (evs(:,3) > 0) .* evs(:,4);
evs(:,6) = evs(:,3) .* evs(:,4);

B = sum(evs(:,5));
%Utilization
U = B / T;
%Average service time
S = U / X;

%Probability of having m parts in the machine (with m = 0, 1, 2)
Y0 = T - B;
Py0 = Y0 / T;

Y1 = min(A(2,1),C(1,1)) - A(1,1);
for m = 2:nC-1
    if A(m,1) >= C(m-1,1)
       Y1 = Y1 + min(A(m+1,1),C(m,1)) - A(m,1);
    end
end
Y1 = Y1 + C(nC,1) - max(A(nA,1),C(nC-1,1));
Py1 = Y1 / T;

Y2 = 0;
if A(2,1)<=C(1,1)
    Y2 = min(A(3,1),C(1,1)) - A(2,1);
end
for m = 3:nC-1
    if A(m,1) <= C(m-1,1) && A(m,1) >= C(m-2,1)
        Y2 = Y2 +min(A(m+1,1),C(m-1,1)) -A(m,1);
    end
end
if A(nA,1) <= C(nC-1,1)
    Y2 = Y2 + C(nC-1,1)-max(C(nC-2,1),A(nA,1));
end
Py2 = Y2 / T;
%Response time less that 30sec & 3min
Pr1 = sum(Rt < 30) / nC;
Pr2 = sum(Rt < 3*60) / nC;

%Inter-arrival time shorter than 1 min
Pi = sum(It < 60) / nA;

%Service time longer than 1 min
si = zeros(nA, 1);
si(1) = AC(1,2) - AC(1,1);
for i = 2 : nA
    si(i) = AC(i,2) - max(AC(i, 1), AC(i-1, 2));
end
Ps = sum(si > 60) / nC;


fprintf(1, "Arrival Rate: %g, Throughput %g\n", Lambda, X);
fprintf(1, "Average inter-arrival time: %g\n", At);
fprintf(1, "Utilization: %g\n", U);
fprintf(1, "Average Service Time: %g\n", S);
fprintf(1, "Average Number of jobs: %g\n", N);
fprintf(1, "Average Response Time: %g\n", R);
fprintf(1, "Probability of having 0/1/2 parts in the machine: %g, %g, %g \n", Py0, Py1,Py2);
fprintf(1, "Probability of having a response time shorter than 30sec and 3min: %g, %g\n", Pr1, Pr2);
fprintf(1, "Probability of having an inter-arrival time shorter than 1 min: %g\n", Pi);
fprintf(1, "Probability of having a service time longer than 1 min: %g\n", Ps);
