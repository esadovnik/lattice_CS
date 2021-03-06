%This function is exactly the same as cosamp, except that it rounds the
%reconstruced signal to the integer lattice at the end of each iteration.
%See cosamp.m for documentation on cosamp

function a = cosampround(phi,u,s,tol,rich,richsteps,signal)
v=u;
[~,N]=size(phi);
a=zeros(N,1);
suppa=uint64.empty;
k=1;

while k < 5                              %Halting criterion option #1
%while norm(v) < tol                      %Halting criterion option #2
%while L^infty norm (phi' * v) < epsilon  %Halting criterion option #3
    
    y = phi'*v;
    omega = uint64(maxindices(y,2*s));
    suppb = union(omega,suppa);
    if rich
        blsq = richardson(phi(:,suppb),u,a(suppb),richsteps);
    else
        blsq = phi(:,suppb)\u;
    end
    b=zeros(N,1);
    b(suppb)=blsq;
    suppa = uint64(maxindices(b,s));
    a = zeros(N,1);
    a(suppa) = b(suppa);
    
    %Here we round the reconstruction to the integer lattice
    a = round(a);
    %a=round(2.0*a)/2.0;
    
    v=u-phi*a;
    
    %residual_round(k)=norm(a-signal);
    k=k+1;
end
%roundingiterations = ['Iterations in rounding = ', num2str(k)];
%disp(roundingiterations)
%plot(residual_round)
end

function maxinds = maxindices(v,s)
Arraycopy = abs(v);
maxinds=zeros(1,s);
for j = 1:s
   [~, maxinds(j)] = max(Arraycopy);
   Arraycopy(maxinds(j)) = 0;
end
end

function z = richardson(tallmat,u,z0,nsteps)
v=tallmat'*u;
z=z0;
for k=1:nsteps
    z=v+z-(tallmat'*(tallmat*z));
end
end