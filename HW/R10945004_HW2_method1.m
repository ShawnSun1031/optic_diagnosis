N = 10000;
ua = 10;
del_z = 0.025;
x = zeros(N,1);
for i = 1:N
    z = 0;
    absorb = false;
    while ~absorb 
        if rand(1) <= ua*del_z
            absorb = true;        
            x(i) = z;
        else
            z = z + del_z;
        end
    end
end
histogram(x,'BinEdge',[0 0:0.025:1 1])
title('Bar graph is the simulated result of r.n.,while the solid line is that from Beer Law');
xlabel('Interval');
ylabel('# of photons absorbs');
        
        
       