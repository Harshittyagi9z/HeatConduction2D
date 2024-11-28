% Plate dimensions
L = 1; % Length of the plate in meters
W = 1; % Width of the plate in meters

% Material properties
alpha = 0.01; % Thermal diffusivity in m^2/s

% Discretization parameters
Nx = 20; % Number of grid points in the x direction
Ny = 20; % Number of grid points in the y direction
dx = L / (Nx - 1); % Grid spacing in x direction
dy = W / (Ny - 1); % Grid spacing in y direction
dt = 0.01; % Time step in seconds

% Simulation time
t_final = 2; % Total simulation time in seconds
n_steps = round(t_final / dt); % Number of time steps

% Initial temperature distribution 
T_init = 20 * ones(Nx, Ny);

% Boundary conditions
T_left = 100; % Left edge temperature (constant)
T_right = 50; % Right edge temperature (constant)
T_top = 25; % Top edge temperature (constant)
T_bottom = 25; % Bottom edge temperature (constant)
% Initialize temperature array
T = T_init;

% Apply boundary conditions to edges (Dirichlet boundary conditions)
T(:, 1) = T_left; % Left boundary
T(:, end) = T_right; % Right boundary
T(1, :) = T_top; % Top boundary
T(end, :) = T_bottom; % Bottom boundary
% Loop for time-stepping (for transient heat conduction)
for t = 1:n_steps
    T_new = T; % Create a copy of the temperature array for the new time step

    % Update interior points (not boundary points)
    for i = 2:Nx-1
        for j = 2:Ny-1
            T_new(i, j) = T(i, j) + alpha * dt * ...
                ((T(i+1, j) - 2*T(i, j) + T(i-1, j)) / dx^2 + ...
                 (T(i, j+1) - 2*T(i, j) + T(i, j-1)) / dy^2);
        end
    end

    % Update the temperature array
    T = T_new;

    % visualize the result at each time step
    if mod(t, 50) == 0 % Display every 50th step
        surf(T);
        title(['Temperature Distribution at Time = ', num2str(t*dt)]);
        xlabel('X Position');
        ylabel('Y Position');
        zlabel('Temperature (°C)');
        colorbar;
        pause(0.01);
    end
end
% Visualizing the final temperature distribution
figure;
surfc(T);
title('Final Temperature Distribution');
xlabel('X Position');
ylabel('Y Position');
zlabel('Temperature (°C)');
colorbar;
