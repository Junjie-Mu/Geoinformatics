clear all;

files = ["Trace1.csv", "Trace2.csv", "Trace3.csv"];
for i = 1:length(files)
    IS = csvread(files{i});
    nA = size(IS, 1);
    %Inter-arrival time
    Ai = IS(:,1);
    %Service time
    Si = IS(:,2);
    %Busy time
    B = sum(Si);
    %Average service time
    S = B / nA;
    %Total time
    T= sum(Ai);
    %Arrival
    A = zeros(nA,1);
    for j = 2:nA
        A(j) = Ai(j - 1) + A(j - 1);
    end
    %Completion
    C = zeros(nA,1);
    C(1,1) = Si(1,1);
    for k = 2:nA
        C(k,1) = Si(k,1) + max(A(k,1) , C(k - 1,1));
    end
    %Average response time
    Rt = C(:,1) - A(:,1);
    W = sum(Rt);
    Lambda = nA / T;
    R = W / nA;
    fprintf(1, "Average response time of %s: %g\n", files(i), R);
    %Utilization
    U = B / T;
    fprintf(1, "Utilization of %s: %g\n", files(i), U);
    %Frequency of idle
    count = 0;
    start = 1;
    for n = 1:nA
        sum1 = sum(Ai(start:n,1));
        sum2 = sum(Si(start:n,1));
        if sum1 > sum2
            count = count + 1;
            start = n+1;
        end
    end
    Fi = count / T;
    fprintf(1, "Frequency of idle of %s: %g\n", files(i), Fi);
    %Average idle time
    Y0 = T - B;
    Ait = Y0 / count;
    fprintf(1, "Average idle time %s: %g\n", files(i), Ait);
    fprintf(1, "=================\n");
end