x_fivetimes = zeros(5,10000);
N = zeros(5,20);
for i=1:5
    x = rand(1,10000);  
    nbins = 20;
    edges = 0:0.05:1;
    Y= discretize(x,edges);
    N(i,:) = histcounts(Y);

end
x_mean = mean(N)
x_std = std(N)
b = histogram(x,nbins);
title('1st run');
xlabel('Interval');
ylabel('Number of random numbers');


