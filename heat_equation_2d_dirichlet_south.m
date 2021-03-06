%Ludovico Fossà 11/2020
%Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)

%% Two dimensional heat equation FD CRANK-NICHOLSON SCHEME
% dirichlet boundary conditions
clc
clear variables
close all

%% INPUT
M=60; %discretization (number of cells per edge)
deltat=1e-2;
ttot=200; %simulation duration
nshow=5;
niter=round(ttot/deltat);
deltasave=round(niter/(nshow-1));
kappa=11.72*1e-4; %thermal diffusivity of iron
%boundary condition - dirichlet (fixed)
%internal initial conditions
Tint_init=273.15;

%% OUTPUT
N=M+1; %discretization (number of nodes)
n=N-2; %inner coordinates
deltax=1/M; %grid step
alpha=kappa*deltat/2/deltax^2; %scheme coefficient
beta1=(1+4*alpha)/alpha;
beta2=(1-4*alpha)/alpha;

%matrices and arrays
A=spdiags(repmat([-1,-1,beta1,-1,-1],(N-2)^2,1),[-N+2,-1,0,1,N-2],(N-2)^2,(N-2)^2);
B=spdiags(repmat([1,1,beta2,1,1],(N-2)^2,1),[-N+2,-1,0,1,N-2],(N-2)^2,(N-2)^2);
T=ones((N-2)^2,1)*Tint_init;
xx=linspace(0,1,N)'; %for the figure

%adapt matrices to boundary conditions
B(n,n+1)=0; %angle n
B((n-1)*n+1,(n-1)*n+1-1)=0; %angle (n-1)*n+1
B(n+1:n:n*(n-2)+1,n+1-1:n:n*(n-1)+1-1)=0; %west
B(2*n:n:n*(n-1),2*n+1:n:n*(n-1)+1)=0; %east
A(n,n+1)=0; %angle n
A((n-1)*n+1,(n-1)*n+1-1)=0; %angle (n-1)*n+1
A(n+1:n:n*(n-2)+1,n+1-1:n:n*(n-1)+1-1)=0; %west
A(2*n:n:n*(n-1),2*n+1:n:n*(n-1)+1)=0; %east

%% BOUNDARY CONDITIONS
% TEMPERATURE DISTRIBUTION AT SOUTH BOUNDARY
T0=273.15; %north-west-east boundary fixed dirichlet
deltaT=100;
%Tdistr=@(x) deltaT*exp(-(x-round((N-1)/2)).^2)+T0; %gaussian
Tdistr=@(x) 100+T0; 
test=ones(N-2,N-2);

for i=1:N-2
    test(i,:)=Tdistr(i);
end
%% INTEGRATION
tt=0;
fprintf('Total number of iterations %d\n',niter)
for k=1:niter
    b=B*T;
    tt=tt+deltat; %time advance
    %boundary conditions: angles
    b(1)=b(1)+2*(Tdistr(2)+T0); %SW
    b(n)=b(n)+2*(Tdistr(N-1)+T0); %SE
    b(n*(n-1)+1)=b(n*(n-1)+1)+4*T0; %NW
    b(n^2)=b(n^2)+4*T0; %NE
    %boundary condition: edges
    b(2:n-1)=b(2:n-1)+2*Tdistr(2:N-3)'; %south
    b(n+1:n:n*(n-2)+1)=b(n+1:n:n*(n-2)+1)+2*T0; %west
    b(2*n:n:n*(n-1))=b(2*n:n:n*(n-1))+2*T0; %east
    b(n*(n-1)+2:n^2-1)=b(n*(n-1)+2:n^2-1)+2*T0; %north
    T=A\b; %solves linear system with the best technique (MATLAB)
    if(mod(k,1000)==0) 
        fprintf('Iteration number k=%d\n',k)
    end
end

XY=reshape(T,n,n);
contourf(XY,'EdgeColor','none')
c=colorbar
c.Label.String = 'Temperature isosurfaces [K]';
c.Label.FontSize = 15;
view(90,-90)