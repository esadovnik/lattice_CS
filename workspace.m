%Set up parameters: d is the number of columns, m is the number of rows, s
%is the sparsity level, enorm is the norm of the noise
d=250;
m=100;
s=30;
enorm=0.5;

%randomly generate the matrix phi, the signal, and the error.
phi = (1/sqrt(m))*random(makedist('Normal'),m,d);
signal=zeros(d,1);
signonzeros= randi([-5,5],s,1);
supp=uint32.empty;
while length(union(supp,uint32.empty))<s
    supp = randi(d,1,s);
end
signal(supp)=signonzeros;
e=random(makedist('Normal'),m,1);
e=enorm*e/norm(e);

%Calculate the observed sample vector, and run cosamp on it
obs=phi*signal+e;
reconstruction=cosamp(phi,obs,s,2*enorm,false,5,signal);
hold off

%print out the norm of the difference between the reconstruction and actual
%signal to the command window
norm(reconstruction-signal)