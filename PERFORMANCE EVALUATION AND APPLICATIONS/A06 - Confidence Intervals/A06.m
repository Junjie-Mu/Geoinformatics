%MU JUNJIE ASSIGNMENT6
M = 1000; 
confidence = 0.95;
confidence_d = norminv ((1+confidence)/2);
relative_error = 0.04;
for i = 1:2
    K = 2;
    Ubar = 0;
    US2 = 0;
    Tbar = 0;
    TS2 = 0;
    A_t=0;
    C_t=0;
    Nbar=0;
    NS2=0;
    
    Rv=0;
    Rv2=0;
    alphaU=1;
    Ui=[];
    Ti=[];
    Ni=[];
    Ri=[];
    RVi=[];
    if(i==1)
        while alphaU > relative_error
            R_t = 0;
            B = 0;
            A_t0 = A_t;
            response_times = zeros(1, M);
            for j = 1: M
                arrival = arr_hyper();
                service = ser_erlang();
                C_t = max(A_t,C_t) + service;
                rt = C_t - A_t;
                response_times(j) = rt;
                R_t = R_t + rt;
                A_t = A_t +arrival;
                B = B + service;
            end
            %Utilization
            total_time = C_t - A_t0;
            Ui(K,1) = B/total_time;
            Ubar = sum(Ui)/K;
            US2 = sum((Ui - Ubar).^2)/(K-1);
            Ci_U = [Ubar - confidence_d * sqrt(US2/K) , Ubar + confidence_d * sqrt(US2/K)];
            alphaU = 2* (Ci_U(1,2) - Ci_U(1,1))/(Ci_U(1,2) + Ci_U(1,1));
            %Throughput
            Ti(K,1) = M /total_time;
            Tbar = sum(Ti)/K;
            TS2= sum((Ti - Tbar).^2)/(K-1);
            Ci_T = [Tbar - confidence_d * sqrt(TS2/K) , Tbar + confidence_d * sqrt(TS2/K)];
            %Averahe numbers of jobs
            Ni(K,1) = R_t /total_time;
            Nbar = sum(Ni)/K;
            NS2= sum((Ni - Nbar).^2)/(K-1);
            Ci_N = [Nbar - confidence_d * sqrt(NS2/K) , Nbar + confidence_d * sqrt(NS2/K)];
            %Average response time
            Ri(K,1) = R_t /M;
            Rbar = sum(Ri)/K;
            RS2= sum((Ri - Rbar).^2)/(K-1);
            Ci_R = [Rbar - confidence_d * sqrt(RS2/K) , Rbar + confidence_d * sqrt(RS2/K)];
            %Variance of the response time
            RVi(K,1) = var(response_times);
            RVbar = sum(RVi)/K;
            RVS2= sum((RVi - RVbar).^2)/(K-1);
            Ci_RV = [RVbar - confidence_d * sqrt(RVS2/K) , RVbar + confidence_d * sqrt(RVS2/K)];
            K = K + 1;
        end
            fprintf(1, "scenario 1\n");
            fprintf(1, " Utilization left bound: %f, Utilization right bound: %f\n",  Ci_U(1,1),Ci_U(1,2));
            fprintf(1, " Throughput left bound: %f, Throughput right bound: %f\n",  Ci_T(1,1),Ci_T(1,2));
            fprintf(1, " Average number of jobs left bound: %f, Average number of jobs right bound: %f\n",  Ci_N(1,1),Ci_N(1,2));
            fprintf(1, " Average response time left bound: %f, Average response time right bound: %f\n",  Ci_R(1,1),Ci_R(1,2));
            fprintf(1, " Vaiance of the response time left bound: %f, Vaiance of the response time right bound: %f\n",  Ci_RV(1,1),Ci_RV(1,2));
            fprintf(1, " Number of batches K to reach the desired accuracy: %i\n",  K);
    end
    if(i==2)
        while alphaU > relative_error
            R_t = 0;
            B = 0;
            A_t0 = A_t;
            response_times = zeros(1, M);
            for j = 1: M
                arrival = arr_exp();
                service = ser_uni();
                C_t = max(A_t,C_t) + service;
                rt = C_t - A_t;
                response_times(j) = rt;
                R_t = R_t + rt;
                A_t = A_t +arrival;
                B = B + service;
            end
            %Utilization
            total_time = C_t - A_t0;
            Ui(K,1) = B/total_time;
            Ubar = sum(Ui)/K;
            US2 = sum((Ui - Ubar).^2)/(K-1);
            Ci_U = [Ubar - confidence_d * sqrt(US2/K) , Ubar + confidence_d * sqrt(US2/K)];
            alphaU = 2* (Ci_U(1,2) - Ci_U(1,1))/(Ci_U(1,2) + Ci_U(1,1));
            %Throughput
            Ti(K,1) = M /total_time;
            Tbar = sum(Ti)/K;
            TS2= sum((Ti - Tbar).^2)/(K-1);
            Ci_T = [Tbar - confidence_d * sqrt(TS2/K) , Tbar + confidence_d * sqrt(TS2/K)];
            %Averahe numbers of jobs
            Ni(K,1) = R_t /total_time;
            Nbar = sum(Ni)/K;
            NS2= sum((Ni - Nbar).^2)/(K-1);
            Ci_N = [Nbar - confidence_d * sqrt(NS2/K) , Nbar + confidence_d * sqrt(NS2/K)];
            %Average response time
            Ri(K,1) = R_t /M;
            Rbar = sum(Ri)/K;
            RS2= sum((Ri - Rbar).^2)/(K-1);
            Ci_R = [Rbar - confidence_d * sqrt(RS2/K) , Rbar + confidence_d * sqrt(RS2/K)];
            %Variance of the response time
            RVi(K,1) = var(response_times);
            RVbar = sum(RVi)/K;
            RVS2= sum((RVi - RVbar).^2)/(K-1);
            Ci_RV = [RVbar - confidence_d * sqrt(RVS2/K) , RVbar + confidence_d * sqrt(RVS2/K)];
            K = K + 1;
        end
            fprintf(1, "scenario 2\n");
            fprintf(1, " Utilization left bound: %f, Utilization right bound: %f\n",  Ci_U(1,1),Ci_U(1,2));
            fprintf(1, " Throughput left bound: %f, Throughput right bound: %f\n",  Ci_T(1,1),Ci_T(1,2));
            fprintf(1, " Average number of jobs left bound: %f, Average number of jobs right bound: %f\n",  Ci_N(1,1),Ci_N(1,2));
            fprintf(1, " Average response time left bound: %f, Average response time right bound: %f\n",  Ci_R(1,1),Ci_R(1,2));
            fprintf(1, " Vaiance of the response time left bound: %f, Vaiance of the response time right bound: %f\n",  Ci_RV(1,1),Ci_RV(1,2));
            fprintf(1, " Number of batches K to reach the desired accuracy: %i\n",  K);

    end
end




%arrival time(hyper-exponential)
function arrival_time = arr_hyper()
    p1 = 0.1;
    lambda1 = 0.02;
    lambda2 = 0.2;
    if rand() < p1
        arrival_time = -log(rand()) /lambda1;
    else
        arrival_time = -log(rand()) /lambda2;
    end
end

%service time(Erlang)
function service_time = ser_erlang()
    lambda_Erlang = 1.5;
    k_Erlang  = 10;
    service_time = sum(-log(rand(1,k_Erlang)) / lambda_Erlang);
end

%arrival time(Exponential)
function arrival_time = arr_exp()
    lambda = 0.1;
    arrival_time = -log(rand())/lambda;
end

%service time(Uniform)
function service_time = ser_uni()
    a = 5;
    b  = 10;
    service_time = a + (b-a)*rand();
end

