clear all
close all
clc

% Alejandro M. Calvo  
% Numerically solves the three-component mixture of sticky and non-sticky EPYC1 and solvent,
% i.e. conservation eqs. of the volume fraction phi_s and phi_ns,
% considering the Flory-Huggins free energy density. We use second-order finite differences to discretize 
% the spatial derivatives and the built-in time-adaptive ode-solve ode45 to integrate the discretized system
% of equations in time. The output in this file is a colormap of sticky EPYC1 volume fraction over time

% Parameters
L      = 25;  % Domain length (square domain)
Nx     = 100; % Number of grid points in x
Ny     = 100; % Number of grid points in y
dx     = L/(Nx-1);
dy     = L/(Ny-1);
dt     = 0.25;  % Time step
tfinal = 150.0; % Final time
vt     = 0:dt:tfinal;
kappa  = 1;   % Gradient energy coefficient
chi    = 6;   % Interaction parameter
k      = 0.1; % Switching rate
D      = 10;  % Dimensionless mobility


% Define my own colormap
map = uint8(flipud([136 8   0
                    174 47  21
                    207 84  50
                    238 136 85
                    246 168 131
                    243 215 185]));

% Set up and writing the movie.
writerObj = VideoWriter('video_coacervation','MPEG-4'); % movie name.
writerObj.FrameRate = 15; % Frames per second. Larger number correlates to smaller movie time duration. 
open(writerObj);


% Initial conditions
phi_s0  = 0.4 + 0.05*rand(Nx,Ny);
phi_ns0 = 0.1 + 0.05*rand(Nx,Ny);

phi0 = [phi_s0(:); phi_ns0(:)];

vx    = 0:dx:L;
vy    = 0:dy:L;
[X,Y] = meshgrid(vx,vy);

% Function for Model B dynamics
model_b = @(t, rho) modelb(rho,Nx,Ny,dx,dy,kappa,chi,k,D);

% Time stepping
for i =1:length(vt)-1

% Solve using ode45
[t,phi] = ode45(model_b,vt(i):dt:vt(i+1),phi0);

% Reshape column vector solution back to the 2D grid
phi_end = phi(end,:)';
phi_s  = reshape(phi_end(1:Nx*Ny), [Nx,Ny]);
phi_ns = reshape(phi_end(Nx*Ny+1:end), [Nx,Ny]);

figure(1);
contourf(X,Y,phi_s,50,'LineColor','none');
axis equal
colormap("winter")
%colormap(map)
colorbar;
set(gca,'FontSize',12)
title('Sticky EPYC1');
xlabel('Coordinate \itx','FontSize',12);
ylabel('Coordinate \ity','FontSize',12);
hcb = colorbar;
hcb.Ruler.Color='w';
clim([0 1])
set(gcf, 'Position',  [400, 400, 800, 400])
%set(gca,'ColorScale','log')
set(gca,'XColor', 'w','YColor','w')
set(gca,'color','k');
set(gcf,'color','k');
drawnow


display(['Time t = ' num2str(vt(i))])

frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
writeVideo(writerObj, frame);


end


close(writerObj);% Saves the movie

% Laplacian discretization
function lap = laplacian(v,dx,dy)
    lap = (circshift(v,[0,-1])- 2*v+circshift(v,[0,1]))/dx^2+...
          (circshift(v,[-1,0])- 2*v+circshift(v,[1,0]))/dy^2;

end

% Free energy minimization (functional derivative), i.e. chemical potential
function [mu_s, mu_ns] = chemical_potentials(phi_s,phi_ns,kappa,chi,dx,dy)

    phi_sol = 1 - phi_s - phi_ns;

    % Avoid log singularities
    eps     = 1e-12;
    phi_s   = max(phi_s, eps);
    phi_ns  = max(phi_ns, eps);
    phi_sol = max(phi_sol, eps);

    % Sticky species (your original)
    mu_s = log(phi_s) - log(phi_sol) + chi*(-2*phi_s) ...
           - kappa*laplacian(phi_s,dx,dy);

    % Passive nonsticky component
    mu_ns = 1 + log(phi_ns)- log(phi_sol);

    % (You may add a coupling term with solvent if desired.)
end


% Spatial discretization of Model B dynamics
function drho_dt = modelb(rho,Nx,Ny,dx,dy,kappa,chi,k,D)

    % extract 2 fields
    phi_s  = reshape(rho(1:Nx*Ny), [Nx,Ny]);
    phi_ns = reshape(rho(Nx*Ny+1:end), [Nx,Ny]);

    % chemical potentials
    [mu_s,mu_ns] = chemical_potentials(phi_s,phi_ns,kappa,chi,dx,dy);

    % model B dynamics (equal mobilities assumed)
    dphi_s  = D*laplacian(mu_s,dx,dy)-k*phi_s + phi_ns;
    dphi_ns = D*laplacian(mu_ns,dx,dy)+k*phi_s - phi_ns;

    % reshape
    drho_dt = [dphi_s(:); dphi_ns(:)];

end
