function T4 = One_D_Trasient_Heat_Conduction_CN_Refractory(fa,fb,Tin,t)
% Diffusion equation - 1D - Crank Nicolson
% Neumann BC conditions - Constant heat flux at boundaries % Unconditionally Stable
% A*T_n+1 = B*T_n+4*d*delta_x *T_BC = C
clear clc

% Inputs
alpha = 2.0e-7; % Thermal diffusivity, m2/s (Heat Conduction - Ozisik & Hahn)
kt = 0.6; % Thermal Conductivity, W/(m K)
%t = 720; % total time, s eg t = 100, 2000
nt = 1000; % total no. of time steps, nt = 1, 200
delta_t=t/nt; % timestep, s
xlength= 0.3; % xlength, m
%ylength = 1; % ylength for graphics prsposes
nx = 20; % total no. of spatial grids slong x direction, eg nx = 5, 20
ny = 1; % total no. of spatial grids along y direction (for graphics purpose)
delta_x = xlength/nx; % delta_x = delta_y, m
%delta_y=ylength/ny; % delta_x = delta_y, m (for graphics purpose)
%fa = -38600; % Heat Flux fa= q = - k* (dT/dx), W/m2 (typ)
%fb = 38600;
%Tin = 100;
% Tmax = max([Tend1,Tend2, Tin]);
% Solution
n = nx + 1; % total no. of nodes/points in the domain
d = alpha*delta_t/delta_x^2; % diffusion number
%x(1) = 0;
for i = 1:n
end
%x;
%x(i) = x(1) + (i-1)* delta_x;
% Creating Initial and Boundary conditions B vector
IC=zeros(n,1);
for k = 1:1
    for i = 1:n
        IC(i,k) = Tin;
    end
end
%IC;
BC=zeros(n,1);
for k = 1:1 
     for i = 1:n
        if i == 1
            BC(i,k) = (fa/kt);
        elseif i == n
            BC(i,k)=(fb/kt);
        else
            BC(i,k) = 0;
        end
    end
end
%BC;
% Creating B matrix
B = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i==j
            B(i,j) = 2*(1 - d);
        elseif (i == j+1) && (i<n)
            B(i,j) = d;
        elseif (i == j+1) && (i == n)
             B(i,j) = 2*d;
        elseif (i == j-1) && (j == 2)
             B(i,j) = 2*d;
        elseif (i == j-1) && (j>2)
            B(i,j) = d;
        else
             B(i,j) = 0;
        end
    end
end
%B;
C=B*IC+ 4*d*delta_x*BC;
% Creating A matrix (see below to save memory by avoiding storing zeros)
%{
A = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i==j
            A(i,j) = 2*(1+d);
        elseif (i == j+1) && (i<n)
            A(i,j) = -d;
        elseif (i == j+1) && (i == n)
            A(i,j) = -2*d;
        elseif (i == j-1) && (i == 1)
            A(i,j) = -2*d;
        elseif (i == j-1) && (i > 1)
            A(i,j) = -d;
            
        else
            A(i,j) = 0;
        end
    end
end
A;
%}

% Algorithm for utilizing Thomas Algorithm
% A matrix diagonal, upper and lower values
dA=zeros(1,n);
uA=zeros(1,n);
lA=zeros(1,n);
for i = 1:n
    dA(i) = 2*(1+d);
    if (i == 1)
        uA(i) = -2*d;
    elseif (i < n) && (i> 1)
        uA(i) = -d;
    else
        uA(i) = 0;
    end
    if (i>1) && (i<n)
        lA(i) = -d;
    elseif (i == n)
        lA(i) = -2*d;
    else
        lA(i) = 0;
    end
end
%lA;
%uA;
%dA;

%{
% creating T vector
T1 = A\C;
%}
T1(1:n,1:ny+1,1)=Tin;
for k = 2:nt+1
    T = Thomas_Algorithm_Func1(dA,lA,uA,C); 
    T=T';
    for j = 1:ny+1
    T1(:,j,k) = T;
    end
    C=B*T+2*d*delta_x*BC;
end

%T1;
%n;
%Tmin = min(T1());
%Tmax = max(T1());
T2 = T1;
% Temperature Profile
%subplot(1,2,1)
% Initial Temperature Profile
%T3 = T2(:,:,1);
%subplot(2,2,1)
%T3 = rot90(T3);
%x= 0:delta_x:xlength;
%y= 0:delta_y:ylength;
%imagesc(x,y,T3)
%colorbar;
%title(['Temperature Profile - ', 'T in deg.C', '@ time(t) = ', num2str(0), 's'])
%xlabel(['x in m'])
%set(gca,'ytick', [])
% Final Temperature Profile
T4 = T2(:,:,nt+1);
%{
subplot(2,2,3)
T4 = rot90(T4);
x = 0:delta_x:xlength;
y= 0:delta_y:ylength;
imagesc(x,y,T4);
colorbar;
title(['Temperature Profile - ', 'T in deg.C', '@ time(t) = ', num2str(t), 's'])
xlabel(['x in m'])
set(gca,'ytick', [])
% Temperature Profile - Animation
T5 = T2(:,:,1)';
[x,y] = meshgrid(0:delta_x:xlength, 0:delta_y-ylength);
x;
y;
for k = 1:nt+1
    subplot(1,2,2)
    colorbar;
   h=surf(x,y,T5,'EdgeColor','none');
    shading interp
    axis ([0 xlength 0 ylength Tmin Tmax]);
    title({['Transient Heat Conduction']:['time (\itt) = ',num2str((k-1)*delta_t), 's']})
    colorbar;
    drawnow;
    pause(0.001);
    refreshdata(h);
    if k~= nt+1
        T5 = T2(:,:,k+1)';
    else
        break;
    end
end
%}
end