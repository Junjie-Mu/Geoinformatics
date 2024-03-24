%MU JUNJIE ASSIGNMENT4
clear all

files = ["Trace1.csv", "Trace2.csv"];
for i = 1:length(files)
    trace = csvread(files{i});
    %1 The first three moments
    Mean = mean(trace);
    Moment2 = mean(trace.^2);
    Moment3 = mean(trace.^3);
    Variance = mean((trace-Mean).^2);
    Deviation = sqrt(Variance);
    Skewness = mean(((trace-Mean)/Deviation).^3);
    %2 Fit Uniform, Exponential, Erlang
    Uniform_a = Mean - 1/2*sqrt(12*(Moment2 - Mean^2));
    Uniform_b = Mean + 1/2*sqrt(12*(Moment2 - Mean^2));

    Exponential_rate = 1/Mean;

    Erlang_stages = round(Mean^2/Variance);
    Erlang_rate = Erlang_stages/Mean;
    %Fit Weibull, Pareto
    wbldist = fitdist(trace, 'Weibull');

    Weibull_shape = wbldist.ParameterValues(2);
    Weibull_scale = wbldist.ParameterValues(1);

    mme_pareto = @(param) [(param(1) * param(2)) / (param(1) - 1) - Mean;
                        (param(1) * param(2)^2) / ((param(1) - 1)^2 * (param(1) - 2)) - Variance];
    param_pareto = fsolve(mme_pareto, [3.0, Mean]);
    Pareto_shape = param_pareto(1);
    Pareto_scale = param_pareto(2);
    %Fit Hyper-Exponential and Hypo-Exponentialdistributions 
    params_hyper = mle(trace, 'pdf', @(X,l1,l2,p1)HyperExp_pdf(X,[l1,l2,p1]), ...
        'Start', [0.8/Mean, 1.2/Mean, 0.4]);
    params_hypo = mle(trace, 'pdf', @(X,l1,l2)HypoExp_pdf(X,[l1,l2]), ...
        'Start', [1/(0.3*Mean), 1/(0.7*Mean)]);

    %Draw empirical CDF of the samples with the CDF
    figure;
    plot(sort(trace),[1:25000]/25000, "-", "DisplayName", "empirical CDF",'LineWidth', 1.5)
    hold on;
    plot([0:7000]/10,HyperExp_cdf([0:7000]/10,[params_hyper(1),params_hyper(2),params_hyper(3)]),"-",'DisplayName', 'Hyper-Exponential','LineWidth', 1.5);
    plot([0:7000]/10,HypoExp_cdf([0:7000]/10,[params_hypo(1),params_hypo(2)]),"-",'DisplayName', 'Hypo-Exponential','LineWidth', 1.5);
    title(['CDF Comparasion of ', files(i)]);
    xlabel('Value');
    ylabel('CDF');
    legend('Location', 'Best');
    hold off;
       
    fprintf("%s:\n", files(i));
    fprintf(' The first three moments are %f, %f, %f\n', Mean, Moment2, Moment3);
    fprintf(' Uniform_a:%f\n Uniform_b:%f\n', Uniform_a, Uniform_b);
    fprintf(' Exponential_rate:%f\n', Exponential_rate);
    fprintf(' Erlang_stages:%f\n Erlang_rate:%f\n', Erlang_stages, Erlang_rate);
    fprintf(' Weibull_shape:%f\n Weibull_scale:%f\n', Weibull_shape, Weibull_scale);
    fprintf(' Pareto_shape:%f\n Pareto_scale:%f\n', Pareto_shape, Pareto_scale);
    fprintf(' Hyper Exponential first rate:%f\n Hyper Exponential second rate:%f\n Hyper Exponential probability of first branch:%f\n', params_hyper(1), params_hyper(2), params_hyper(3));
    fprintf(' Hypo Exponential first rate:%f\n Hypo Exponential second rate:%f\n', params_hypo(1), params_hypo(2));

end