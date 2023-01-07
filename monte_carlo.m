num_monte_carlo = 100;
RMSE = zeros(num_monte_carlo,1);
model.plot = false;
for i=1:1:num_monte_carlo
    main;
    RMSE(i) = sqrt(1/length(truth)*sum(diag((truth(1:2,:) - est.x(1:2,:))*(truth(1:2,:) - est.x(1:2,:))')));
end
num_monte_carlo = 100;
bar_RMSE = sum(RMSE)/num_monte_carlo;