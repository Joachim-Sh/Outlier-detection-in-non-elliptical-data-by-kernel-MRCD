function data = halfkernel(N1,N2, minx, r1, r2, noise, ratio)


phi1 = rand(N1,1) * pi;
inner = [minx + r1 * sin(phi1) - .5 * noise  + noise * rand(N1,1) r1 * ratio * cos(phi1) - .5 * noise + noise * rand(N1,1) ones(N1,1)];
    
phi2 = rand(N2,1) * pi;
outer = [minx + r2 * sin(phi2) - .5 * noise  + noise * rand(N2,1) r2 * ratio * cos(phi2) - .5 * noise  + noise * rand(N2,1) zeros(N2,1)];

data = [inner; outer];
end