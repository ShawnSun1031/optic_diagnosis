%% Bar graph of simulate result
N = 10000;
ua = 10;
del_z = 0.025;
for epoch = 1:5 
    x = zeros(N,1);
    for i = 1:N
        s = -log(rand(1))/ua;
        x(i) = floor(s/del_z)*del_z;
    end
subplot(3,2,epoch);histogram(x,'BinEdge',[0 0:0.025:1 1])
title('Bar graph is the simulated result of r.n.,while the solid line is that from Beer Law');
xlabel('Interval');
ylabel('# of photons absorbs');
%% Beer Law
x_beer = [0:0.025:1];
I0 = N;
y = zeros(1,41);
for i = 1:41
    abs = I0 - I0*exp(-ua*del_z);
    y(i) = abs;
    I0 = I0*exp(-ua*del_z);
end
hold on
plot(x_beer,y)

end



