%This function is exactly the same as cosamp, except that it rounds the
%reconstruced signal to the integer lattice only at the end of the whole process.
%See cosamp.m for documentation on cosamp

function aa = cosamproundend(phi,u,s,tol,rich,richsteps,signal)
v=u;
[~,N]=size(phi);
aa=zeros(N,1);
suppa=uint64.empty;
k=1;

while k < 5                        %Halting criterion option #1
%while norm(v) < tol                      %Halting criterion option #2
%while L^infty norm (phi' * v) < epsilon  %Halting criterion option #3

    y = phi'*v;
    omega = uint64(maxindices(y,2*s));
    suppb = union(omega,suppa);
    if rich
        blsq = richardson(phi(:,suppb),u,aa(suppb),richsteps);
    else
        blsq = phi(:,suppb)\u;
    end
    b=zeros(N,1);
    b(suppb)=blsq;
    suppa = uint64(maxindices(b,s));
    aa = zeros(N,1);
    aa(suppa) = b(suppa);
    v=u-phi*aa;
    
    %residual(k)=norm(a-signal);
    k=k+1;
end

aa = round(aa);
%aa = round(aa*2.0)/2.0;

%residual(end)=norm(a-signal);
%endroundingiterations = ['Iterations in rounding at end= ', num2str(k)];
%disp(endroundingiterations)
%plot(residual)
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