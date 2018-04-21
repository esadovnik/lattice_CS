%cosamp is called with 7 arguments:
%
%  phi: The sampling matrix
%
%  u: the observed samples
%
%  s: the sparsity level
%
%  tol: the tolerance level. If the approximated signal ever would give a
%               sample vector whose difference with the observed sample
%               vector has norm smaller than tol, the iteration stops and
%               the current approximated signal is returned as the
%               reconstructed signal
%
%  rich: a true/false value -- true corresponds to using richsteps steps of
%               richardson iteration to find the least squares solution at
%               each step, false corresponds to using matlab's built-in
%               least squares solver. Have not found a regime in which
%               richardson iteration does better (and it diverges in
%               some regimes, whereas matlab's built in function never
%               does)
%
%  richsteps: the number of steps to be used for richardson iteration. If
%               rich is false, then this is not used in the algorithm
%
%  signal: the true signal. This is only used to generate a plot of the
%               residual at each step, not in the reconstruction algorithm
%               itself

function a = cosamp(phi,u,s,tol,rich,richsteps,signal)
%  initialize v, a, N, suppa, and k

%  v is the unreconstructed part of the observed samples at each step
v=u;

%  N is the number of columns of phi, larger than the number of rows
[~,N]=size(phi);

%  a is the current reconstructed signal. It is initialized to zero.
a=zeros(N,1);

%  suppa is the support of a, which we need to keep track of. Since a is
%               initialized to zero, suppa must be initialized to be empty
suppa=uint64.empty;

%  k is the iteration number.
k=1;

while norm(v)>tol && k < 1000
    % first we form the signal proxy y and the indices of its 2s largest
    % entries
    y = phi'*v;
    omega = uint64(maxindices(y,2*s));
    
    % Then we merge the support of the 2s largest entries with the support
    % of the reconstructed signal from the previous step, and perform least
    % squares to get the best possible reconstruction b supported on the
    % merged set
    suppb = union(omega,suppa);
    if rich
        blsq = richardson(phi(:,suppb),u,a(suppb),richsteps);
    else
        blsq = phi(:,suppb)\u;
    end
    b=zeros(N,1);
    b(suppb)=blsq;
    
    % We threshold b to its s largest entries to get the current signal
    % reconstruction, and update the support of the reconstruction.
    suppa = uint64(maxindices(b,s));
    a = zeros(N,1);
    a(suppa) = b(suppa);
    
    % Finally we update the current samples
    v=u-phi*a;
    
    % We record the residual from the actual signal for plotting later, and
    % increment the iteration counter
    residual(k)=norm(a-signal);
    k=k+1;
end

%  We display the number of iterations performed and the plot of the
%  residual vs iteration count, and return the last reconstructed signal
disp(k)
plot(residual)
end

%  The function maxindices takes in a vector v and an integer s and finds
%  the indices of the s largest (in absolute value) components of v. Note:
%  the returned set of indices is listed in the order of decreasing
%  absolute magnitude of the corresponding entries of v, NOT in ascending
%  order.

function maxinds = maxindices(v,s)
Arraycopy = abs(v);
maxinds=zeros(1,s);
for j = 1:s
   [~, maxinds(j)] = max(Arraycopy);
   Arraycopy(maxinds(j)) = 0;
end
end

% This is the richardson iteration function. It takes in a tall matrix
% tallmat, vectors u and z0, and an integer nsteps. It uses nsteps
% iterations to solve the system tallmat*z=u for z, and uses the vector z0
% as the initial guess.

function z = richardson(tallmat,u,z0,nsteps)
v=tallmat'*u;
z=z0;
for k=1:nsteps
    z=v+z-(tallmat'*(tallmat*z));
end
end