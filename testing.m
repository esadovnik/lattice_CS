%Parameters
d = 250;
m = 100;
s = 10;

%Initialize lists for plotting at the end
enorm_list = [0:0.1:5];
normal_data = zeros(1,51);
end_rounding_data = zeros(1,51);
rounding_data = zeros(1,51);

for enorm = 0:0.1:5
    
    %Initialize lists for averaging
    normal_results = zeros(1,50);
    end_rounding_results = zeros(1,50);
    rounding_results=  zeros(1,50);
    tol = max(enorm, 0.1);
    
    for j=1:50
        
        %Randomly generate the matrix phi, the signal, and the error.
        phi = (1/sqrt(m))*random(makedist('Normal'),m,d);
        signal=zeros(d,1);
        signonzeros= randi([-1,1],s,1);
        supp=uint32.empty;
        while length(union(supp,uint32.empty))<s
            supp = randi(d,1,s);
        end
        signal(supp)=signonzeros;
        e=random(makedist('Normal'),m,1);
        e=enorm*e/norm(e);
        
        %Make the sample and run the three CoSaMPs to get the reconstructions
        u = phi * signal + e;
        normal_reconstruction = cosamp(phi, u, s, tol, false, 5, signal);
        end_rounding_reconstruction = cosamproundend(phi, u, s, tol, false, 5, signal);
        rounding_reconstruction = cosampround(phi, u, s, tol, false, 5, signal);
        
        %Add the reconstructions to the lists for averaging
        normal_results(j) = norm(normal_reconstruction - signal);
        end_rounding_results(j) = norm(end_rounding_reconstruction - signal);
        rounding_results(j) = norm(rounding_reconstruction - signal);
        
    end
        
    %Add average of the 50 trials into the list for plotting
    normal_data(round(10*enorm + 1)) = mean(normal_results);
    end_rounding_data(round(10*enorm + 1)) = mean(end_rounding_results);
    rounding_data(round(10*enorm + 1)) = mean(rounding_results);
    
end

%Plot everything
hold on
plot(enorm_list, normal_data, 'LineWidth',2)
plot(enorm_list, end_rounding_data, 'LineWidth', 6)
plot(enorm_list, rounding_data,'LineWidth', 2)
legend('Normal CoSaMP','Rounding at end','Intermediate rounding')
xlabel('Noise norm')
ylabel('Average error in reconstruction (50 trials)')
        