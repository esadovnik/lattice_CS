%Parameters
d = 250;
m = 100;
s = 10;
enorm_min = 0;
enorm_max = 5;
enorm_step = 0.1;
trials = 50;

%Initialize lists for plotting at the end
enorm_list = enorm_min:enorm_step:enorm_max;
normal_data = zeros(1,round((enorm_max-enorm_min)/enorm_step + 1));
end_rounding_data = zeros(1, round((enorm_max-enorm_min)/enorm_step + 1));
rounding_data = zeros(1, round((enorm_max-enorm_min)/enorm_step + 1));

for enorm = enorm_min:enorm_step:enorm_max
    
    %Initialize lists for averaging
    normal_results = zeros(1,trials);
    end_rounding_results = zeros(1,trials);
    rounding_results=  zeros(1,trials);
    tol = max(enorm, 0.1);
    
    for j=1:trials
        
        %Randomly generate the matrix phi, the signal, and the error.
        phi = (1/sqrt(m))*random(makedist('Normal'),m,d);
        signal=zeros(d,1);
        signonzeros= randi([-1,1],s,1);
        supp=uint32.empty;
        while length(union(supp,uint32.empty))<s
            supp = randi(d,1,s);
        end
        signal(supp)=signonzeros;
        e=rand(m,1);
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
        
    %Add average of the trials into the list for plotting
    normal_data(round(enorm_list,2) == round(enorm,2)) = mean(normal_results);
    end_rounding_data(round(enorm_list,2) == round(enorm,2)) = mean(end_rounding_results);
    rounding_data(round(enorm_list,2) == round(enorm,2)) = mean(rounding_results);
    
end

%Plot everything
plot(enorm_list, normal_data, 'LineWidth',2)
hold on
plot(enorm_list, end_rounding_data, 'LineWidth', 2)
plot(enorm_list, rounding_data,'LineWidth', 2)
hold off
legend('Normal CoSaMP','Rounding at end','Intermediate rounding')
xlabel('Noise norm')
ylabel(['Average error in reconstruction (' num2str(trials) ' trials)'])
        