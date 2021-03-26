% N-S equation using vorticity transport equation


clc 
clear all
%%initialize variable and domain
L_x=0.1;
L_y=0.1;
delta_t=0.001;
 
U=1;
mu=0.001;
tol_phi=10^(-6);   %//error threshold for phi equations -6
tol_w=10^(-2);     %//error threshold for vorticity equation -5
 
N_x=51;
delta_x= L_x/(N_x-1);
N_y=N_x;
delta_y=delta_x;
 
a_x=(mu*delta_t/(delta_x)^2);
a_y=(mu*delta_t/(delta_y)^2);
b_x=delta_t/delta_x;
b_y=delta_t/delta_y;
 
 
 
 
%% 
 
%%initial & boundary conditions for variables
%boundary points
u(1:N_x,N_y,1)=1;
u(1:N_x,1,1)=0;
u(1,1:N_y-1,1)=0;
u(N_x,1:N_y-1,1)=0;
v(1:N_x,N_y,1)=0;
v(1:N_x,1,1)=0;
v(1,1:N_y-1,1)=0;
v(N_x,1:N_y-1,1)=0;
% u(1:N_x,N_y,1)=0;
% u(1:N_x,1,1)=0;
% u(1,1:N_y-1,1)=1;
% u(N_x,1:N_y-1,1)=0;
% v(1:N_x,N_y,1)=0;
% v(1:N_x,1,1)=0;
% v(1,1:N_y-1,1)=0;
% v(N_x,1:N_y-1,1)=0;
%interior points
 
 
u(2:N_x-1,2:N_y-1,1)=0;
v(2:N_x-1,2:N_y-1,1)=0;
 
 
% initial condition for phi
% boundary points
phi(1:N_x,N_y,1)=0;
phi(1:N_x,1,1)=0;
phi(1,1:N_y-1,1)=0;
phi(N_x,1:N_y-1,1)=0;
% interior points
 
phi(2:N_x-1,2:N_y-1,1)=0;
 
 
% initialize w
% boundary points
w(1,1:N_y,1) = -1*phi(2,1:N_y,1)*(2/(delta_x)^2);
w(N_x,1:N_y,1) = -1*phi(N_x-1,1:N_y,1)*(2/(delta_x)^2);
 
w(1:N_x,1,1) = -1*phi(1:N_x,2,1)*(2/(delta_y)^2);
w(1:N_x,N_y,1) = -1* phi(1:N_x,N_y-1,1)*(2/(delta_y)^2)-(2*U/delta_y);
%interior points
w(2:N_x-1,2:N_y-1,1)=0;
err_w=1;
n=1;
w(:,:,2)=w(:,:,1);
iter_w=1;
%% 

while err_w > tol_w
    n=n+1;
    iter_w=iter_w+1;
    
    for j=2:N_y-1
        for i=2:N_x-1
              
            % nodal equation for interior points
              w(i,j,n)=((1-2*(a_x+a_y))*w(i,j,n-1))+((-u(i,j,n-1)*(b_x/2)+a_x)*w(i+1,j,n-1))+((u(i,j,n-1)*(b_x/2)+a_x)*w(i-1,j,n))+((-v(i,j,n-1)*(b_y/2)+a_y)*w(i,j+1,n-1))+((v(i,j,n-1)*(b_y/2)+a_y)*w(i,j-1,n));
            
        end
    end
%% 
% finding phi using poisson equation
 
 
err_phi=1;
m=.9;
iter_phi=1;
phi_n(:,:,n-1)=phi(:,:,n-1);
%for k=1:500
while err_phi>tol_phi
    iter_phi=iter_phi+1;
    for j=2:N_y-1
        for i=2:N_x-1
            
            
            
           
            phi_n(i,j,n-1)=phi(i,j,n-1) + (m/(2*(delta_x^2+delta_y^2)))*((w(i,j,n)*((delta_x*delta_y)^2))+((delta_y^2)*(phi(i+1,j,n-1)+phi(i-1,j,n-1)))+((delta_x^2)*(phi(i,j+1,n-1)+phi(i,j-1,n-1)))+(-2*((delta_x)^2+(delta_y)^2)*phi(i,j,n-1)));
            
        end 
        
            
    end
    
    
    err_phi=abs(max(max((phi_n(:,:,n-1)-phi(:,:,n-1)))));
    phi_n(1:N_x,N_y,n-1)=0;
    phi_n(1:N_x,1,n-1)=0;
    phi_n(1,1:N_y-1,1)=0;
    phi_n(N_x,1:N_y-1,1)=0;
    

phi(:,:,n-1)=phi_n(:,:,n-1);
phi(:,:,n)=phi(:,:,n-1);
end
%% finding u and v from phi
 
for j=2:N_y-1
    for i=2:N_x-1
        u(1:N_x,N_y,n)=1;
        u(1:N_x,1,n)=0;
        u(1,1:N_y-1,n)=0;
        u(N_x,1:N_y-1,n)=0;
        v(1:N_x,N_y,n)=0;
        v(1:N_x,1,n)=0;
        v(1,1:N_y,n)=0;
        v(N_x,1:N_y,n)=0;
        u(i,j,n)=(phi(i,j+1,n)-phi(i,j-1,n))/(2*delta_y);
        v(i,j,n)=-(phi(i+1,j,n)-phi(i-1,j,n))/(2*delta_x);
    end
end
    
              
              w(1,1:N_y,n) = -phi(2,1:N_y,n)*(2/(delta_x)^2);
              w(N_x,1:N_y,n) = -phi(N_x-1,1:N_y,n)*(2/(delta_x)^2);
              w(1:N_x,1,n) = -phi(1:N_x,2,n)*(2/(delta_y)^2);
              w(1:N_x,N_y,n) =- phi(1:N_x,N_y-1,n)*(2/(delta_y)^2)-(2*U/delta_y);
              w(:,:,n+1)=w(:,:,n);
              err_w(n)=(abs(max(max((w(1:N_x,1:N_y,n)-w(1:N_x,1:N_y,n-1))))));
 
end
 
 


%% ploting of velocities
Y=[1:51];
for i=1:max(size(Y))
    y(i)=((Y(i)-1)*delta_y);
end
figure(1)
plot(u(25,:,n),y(:))
            
X=[1:51];
for i=1:max(size(X))
    x(i)=((X(i)-1)*delta_x);
end
figure(2)
plot(x(:),v(:,25,n))

%%contour of x-dir velocity
[X,Y] = meshgrid(x,y);
contourf(X,Y,flip(rot90(rot90(rot90(u(:,:,235)))),2),50,'ShowText','off')
xlabel('x');
ylabel('y');
title('X-velocity Contour') 
colormap("jet");
colorbar
%%
%contour of velocity vector
figure(3)
quiver(X,Y,u(:,:,235),v(:,:,235),3,'r');
xlabel('x');
ylabel('y');
title('Velocity Vector') 




