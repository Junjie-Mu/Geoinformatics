%MU JUNJIE ASSIGNMENT3
clear all

files = ["Trace1.csv", "Trace2.csv", "Trace3.csv"];
for i = 1:length(files)
    trace = csvread(files{i});
    %Mean and the second,third and fourth moments
    Mean = mean(trace);
    moment_2 = mean(trace.^2);
    moment_3 = mean(trace.^3);
    moment_4 = mean(trace.^4);
    %Variance and the third, fourth centered moments
    variance = mean((trace-Mean).^2);
    centered_mom_3 = mean((trace-Mean).^3);
    centered_mom_4 = mean((trace-Mean).^4);
    %Skewness and forth standardized moment
    %Standard Deviation, Coefficient of Variation and Excess Kurtosis
    deviation = sqrt(variance);
    skewness = mean(((trace-Mean)/deviation).^3);
    standard_mom_4 = mean(((trace-Mean)/deviation).^4);

    coefficient = variance / Mean;
    excess_kurtosis = standard_mom_4 - 3;
    %Median, the first and third quartile
    median_value = median(trace);
    quartile_1 = quantile(trace, 0.25);
    quartile_3 = quantile(trace, 0.75);

    fprintf("%s:\n", files(i));
    fprintf(' Mean: %f\n Second moment: %f\n Third moment: %f\n Fourth moment: %f\n', ...
        Mean, moment_2, moment_3, moment_4 );
    fprintf(' Variance: %f\n Third centered moment: %f\n Fourth centered moment: %f\n', ...
        variance, centered_mom_3, centered_mom_4);
    fprintf(' Skewness: %f\n Fourth standardized moment: %f\n', ...
        skewness,  standard_mom_4);
    fprintf(' Standard Deviation: %f\n Coefficient of Variation: %f\n Excess Kurtosis: %f\n', ...
        deviation,  coefficient, excess_kurtosis);
    fprintf(' Median: %f\n First_quartile: %f\n Third quartile: %f\n', ...
        median_value,  quartile_1, quartile_3);

    %Draw the Pearson Correlation Coedfficient
    lags = 1:100;
    correlation_coefficients = zeros(1, 100);
    for m = lags
        X = trace(1:end-m); 
        Y = trace(m+1:end);  
        covariance = mean(X .* Y) - mean(X) * mean(Y);
        deviation_X = std(X);
        deviation_Y = std(Y);
        correlation_coefficients(m) = covariance / (deviation_X * deviation_Y);
    end
    figure;
    plot(lags, correlation_coefficients);
    title(['Pearson Correlation Coefficient of' files(i)]);
    xlabel('Lag (m)');
    ylabel('Correlation Coefficient');

    %Draw the approximated CDF
    sorted_trace = sort(trace);
    n = length(trace);
    row = ((1:n)')/n;
    figure;
    plot(sorted_trace,row);
    title(['Approximated CDF of' files(i)]);
end


