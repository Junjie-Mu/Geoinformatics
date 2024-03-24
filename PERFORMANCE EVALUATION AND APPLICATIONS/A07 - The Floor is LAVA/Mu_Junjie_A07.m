clear all;
%MU JUNJIE ASSIGNMENT7
%Winning probability
pro_win = 0.7*0.8*(0.5*0.25+0.5*0.6)*0.6+0.3*0.3*0.6;
%Average duration
s0 = 1;
s = s0;
t = 0;
T = 10000;
TS1 = 0;
TS2 = 0;
TS3 = 0;
C = 0;
v = 0;

while t < T
    %Entrance
	if s == 1 
        p1 = rand();
        if p1 < 0.7
            p2 = rand();
            if p2 < 0.8
                ns = 2;
                dt = getErlang(4,1.5);
            else
                ns = 4;
                dt = getExp(0.5);
            end
        else
            p3 = rand();
            if p3 < 0.7
                ns = 4;
                dt = getExp(0.25);
            else
                ns = 3;
                dt = getUniform(3,6);
            end
        end
        TS1 = TS1 + dt;
	end

    %C1
	if s == 2 
        p4 = rand();
        if p4 < 0.5
        p5 = rand();
            if p5 < 0.75
                ns = 4;
                dt = getExp(0.4);
            else
                ns = 3;
                dt = getErlang(3,2);
            end
        else
        p6 = rand();
            if p6 < 0.4 
                ns = 4;
                dt = getExp(0.2);
            else
                ns = 3;
                dt = getExp(0.15);
            end 
        end
        TS2 = TS2 + dt;
    end
        
    %C2
    if s == 3
        p7 = rand();
        if p7 < 0.4
            ns = 4;
            dt = getErlang(5,4);
        else
            ns = 5;
            dt = getErlang(5,4);
        end
        TS3 = TS3 + dt;
    end

    %Fall
    if s == 4
        ns = 1;
        dt = 5;
        C = C + 1;
    end

    %Exit
    if s == 5
        ns = 1;
        dt = 5;
        C = C + 1;
    end
	s = ns;
	t = t + dt;
end

duration_avg = (TS1 + TS2 + TS3) / C;
X = C / t * 60;

fprintf(1, "Winning probability: %f\n",  pro_win);
fprintf(1, "Average duration: %f\n",  duration_avg);
fprintf(1, "Throughput(per hour): %f\n",  X);
    

%Erlang
function erlang = getErlang(k_Erlang, lambda_Erlang)
    erlang = sum(-log(rand(1,k_Erlang)) / lambda_Erlang);
end

%Exponential
function lambda = getExp(lambda)
    lambda = -log(rand())/lambda;
end

%Uniform
function uniform = getUniform(a,b)
    uniform = a + (b-a)*rand();
end