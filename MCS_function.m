
%        Monte Carlo simulation of light transport through tissue         

function [coordinate_store, path_store, Rdr, Rd, delta_r, idx, edges,...
    Td, R_unscat, T_unscat, Tdr, R_layers, T_layers, Escaped_bounds,...
    Roulette_weight, R_e, T_e, all_paths, abs_coords, abs_weight] = MCS_function(simulation, tissue, boundaries)

% Make some loaded boundaries global for passing to the photon class
global max_steps, max_steps = boundaries.max_steps;
global max_events, max_events = boundaries.max_events;
global max_radius, max_radius = boundaries.max_radius;
global max_depth, max_depth = sum(tissue.layers(1:end));
global max_layer_number, max_layer_number = length(tissue.layers);

% ----- Set more parameters -----
tissue.muT = tissue.muA + tissue.muSr;
clear_layer = 0; n1 = 1; % Refractive index of air %1.00027717; 
n2 = tissue.refractive_index(1); n3 = 1; % Set n3 for ambient medium

% ----- Storage allocation -----
path_store = []; coordinate_store = {}; R_unscat = 0; T_unscat = 0;
all_paths = []; abs_coords = {}; abs_weight = {};
% Radial diffuse reflectance % radial bins spacing 0.09 mm
Nr = ceil(max_radius/0.09); delta_r = max_radius/Nr; % Capture radius/Bin number 
Rdr = zeros(2,Nr); idx = zeros(1,Nr);
for i = 0:Nr-1   % Wang
    idx(1,i+1) = ((i+0.5)+(1/(12*(i+0.5))))*delta_r; % indices
end
for i = 0:Nr
    r = ((i)+(1/(12*(i))))*delta_r; % edges
    Rdr(1,i+1) = r;
end
Rdr(1,1) = 0; % when using edges
Tdr = Rdr;
edges = Rdr(1,:); % Rdr(1,1:end-1) = idx;

% PSF of diffuse reflectance
N = ceil(max_radius/0.09); b = linspace(0,max_radius,N); 
b = [fliplr(-b),b(2:end)]; Rd = zeros(length(b),length(b),2);
for i = 1:length(b)
   gx = b(i);
   for j = 1:length(b)
       gy = b(j); gr = sqrt(gx^2 + gy^2); Rd(i,j,1) = gr;
   end
end
Td = Rd;

layer_edges = zeros(1,length(tissue.layers)+1);
for n = 1:length(tissue.layers)
    layer_edges(n+1) = sum(tissue.layers(1:n));
end

R_e = struct(); R_e.c = [0 0 0]; R_e.p = [0]; T_e = R_e;

R_layers = zeros(1,size(tissue.layers,2)); % Internal reflectance in layer
T_layers = R_layers; % Transmittance from each layer
Escaped_bounds = {}; % Exceeded max radius
Escaped_bounds.coordinates = [0 0 0];
Escaped_bounds.weight = [0];
Escaped_bounds.r = [0];
Roulette_weight = 0; % Weight added to the n*1 photon weight

% ----- Run the procedure -----
rng('shuffle'); % Shuffle random numbers
%h = waitbar(0,'Running simulation...');
for n = 1:simulation.number_of_photons

    % ----- Initialise a new photon -----
    p = photon(); p.number = n;
    photon_store = [p.coordinates]; % Log initial position
    p.path = zeros(1,length(tissue.layers)); % Record path in each layer    
    
    % ----- Specular reflectance -----
    if clear_layer == 0 
        Rsp = (n1-n2)^2/(n1+n2)^2;                
    else % If the photon travels through a clear medium   
        r1 = (n1-n2)^2/(n1+n2)^2; r2 = (n3-n2)^2/(n3+n2)^2;
        Rsp = r1 + ((1-r1)^2*r2)/(1-r1*r2);
    end            
    
    while(p.active == 1) 
        if p.step_number == 0 % Specular reflectance decreases the weight      
            p.weight = p.weight - Rsp; 
        end        
        if p.step_size == 0 % Reset step size for new or scattered photons
            p.step_size = -log(rand);        
        end     
        
        % ----- Layer boundaries -----
        % A photon can be internally reflected or transmit across boundary
        I = layer_edges(p.layer_number); J = layer_edges(p.layer_number + 1); 
        
        % Distance from photon location to the current layer boundaries I,J
        if p.direction(3) < 0
            distance = (I-p.coordinates(3))/p.direction(3);
        elseif p.direction(3) == 0
            distance = inf;
        elseif p.direction(3) > 0
            distance = (J-p.coordinates(3))/p.direction(3);
        end

        % Check if the dimensionless stepsize is greater than the distance
        interaction = distance*tissue.muT(p.layer_number) <= p.step_size;
        if interaction % If true, the photon hits the tissue boundary 
            % Photon is moved to the boundary and the step size updated
            p.coordinates = p.coordinates + p.direction*distance;                                                                            
            p.step_size = p.step_size - distance*tissue.muT(p.layer_number);
            p.path(p.layer_number) = p.path(p.layer_number) + distance;                                 
            if p.active ~= 0 % If still within boundaries

            % ----- Internal reflection vs transmittance -----
            alpha_i = acosd(abs(p.direction(3)));   
            nI = tissue.refractive_index(p.layer_number); % Current layer
            if (p.layer_number == 1 && sign(p.direction(3)) == -1) ||...
            (p.layer_number == length(tissue.layers) && sign(p.direction(3)) == 1)
                nT = 1; % Refractive index of air               
            else 
                nT = tissue.refractive_index(p.layer_number+sign(p.direction(3)));
            end            
            alpha_t = asind((nI*sind(alpha_i))/nT); % Snell's law 
            if alpha_i > asind(nT/nI)
                internal_reflectance = 1;
            else % Fresnel's formulas:
                internal_reflectance = (1/2)*((sind(alpha_i-alpha_t))^2/...
                ((sind(alpha_i+alpha_t))^2)+((tand(alpha_i-alpha_t))^2)/...
                ((tand(alpha_i+alpha_t))^2));
            end          
            
            % Decide if the photon is internally reflected or transmitted
            if rand <= internal_reflectance % Internally reflected
                p.direction = p.direction.*[1 1 -1]; % Reverse z-direction
                R_layers(p.layer_number) = R_layers(p.layer_number) + p.weight;
            else % Transmitted:
                T_layers(p.layer_number) = T_layers(p.layer_number) + p.weight;
                dir = sign(p.direction(3)); % +/- is in/out wrt surface 
                p.layer_number = p.layer_number + dir; 
                if p.active == 2 % Photon escaped the tissue
                    if p.layer_number == 0 
                        p.coordinates(3) = 0; % Already 0, rounding
                    elseif p.layer_number == (length(tissue.layers)+1)
                        p.coordinates(3) = layer_edges(end);
                    end
                    nT = 1; % Refractive index of air
                else 
                    nT = tissue.refractive_index(p.layer_number); % To 
                end 
                % Find the new directions
                nI = tissue.refractive_index(p.layer_number - dir); % From
                x_dir = p.direction(1)*nI/nT;
                y_dir = p.direction(2)*nI/nT;
                z_dir = sign(p.direction(3))*cosd(alpha_t);   
                p.direction = [x_dir, y_dir, z_dir]; 
            end 
            end
        else % If not hitting the layer boundary, the photon moves by s/muT
            % ----- Photon movement -----            
            p.coordinates = p.coordinates + (p.direction)*p.step_size/...
                tissue.muT(p.layer_number);
            p.path(p.layer_number) = p.path(p.layer_number) + ...
                p.step_size/tissue.muT(p.layer_number); 
            p.step_number = p.step_number + 1; p.step_size = 0; % Reset
            if p.active ~= 0 % If still within boundaries

            % ----- Photon absorption -----      
            weight_change = (tissue.muA(p.layer_number)/...
                tissue.muT(p.layer_number))*p.weight;
            p.weight = p.weight - weight_change;
	    abs_coords = [abs_coords;p.coordinates];
	    abs_weight = [abs_weight;weight_change];                                               

            % ----- Photon scattering -----                                    
            if tissue.g == 0 % Isotropic scattering case
                cos_theta = 2*rand -1;
            else % Anisotropic case, Henyey-Greenstein phase function
                g = tissue.g(p.layer_number);
                cos_theta = (1/(2*g))*(1+g^2-((1-g^2)/(1-g+2*g*rand))^2);
            end
            theta = acos(cos_theta); % The deflection angle
            phi = 2*pi*rand; % The azimuthal angle is uniformly distributed 
            % Calculate the new direction based on these angles: [1]
            if abs(p.direction(3)) < 0.99999 
                x_dir = sin(theta)*(p.direction(1)*p.direction(3)*cos(phi)-...
                    p.direction(2)*sin(phi))/(sqrt(1-p.direction(3)^2))+...
                    p.direction(1)*cos(theta);
                y_dir = sin(theta)*(p.direction(2)*p.direction(3)*cos(phi)+...
                    p.direction(1)*sin(phi))/(sqrt(1-p.direction(3)^2))+...
                    p.direction(2)*cos(theta);
                z_dir = -sin(theta)*cos(phi)*sqrt(1-p.direction(3)^2)+...
                    p.direction(3)*cos(theta);
            else % Photon direction is close to the z-axis [Wang 3.25] 
                x_dir = sin(theta)*cos(phi);                                
                y_dir = sin(theta)*sin(phi);
                z_dir = sign(p.direction(3))*cos(theta);
            end
            p.direction =  [x_dir y_dir z_dir];
            end
        end
        if p.active == 1 % Russian roulette
            if p.weight < boundaries.threshold_weight % Small weight 
                if rand <= 1/boundaries.m
                    Roulette_weight = Roulette_weight + (boundaries.m*p.weight - p.weight);
                    p.weight = boundaries.m*p.weight; 
                else
                    Roulette_weight = Roulette_weight + p.weight;
                    p.weight = 0;
                end
            end
        end            
        % Store the coordinates for each photon step
        photon_store = [photon_store ; p.coordinates];              
    end % Photon no longer propagating
    coordinate_store{end+1} = photon_store;
    
    % --- Update with results from the last photon ---
    all_paths = [all_paths ; p.path]; % From all photons
    if p.active == 2 % Photon escaped the tissue        
        % Store the path lengths of remitted photons
        if p.layer_number == 0
            path_store = [path_store ; p.path]; 
        end	
        % Score escaped photon weight into radial reflectance array
        R = max_radius*max_radius + 1; x = -1;
        r = sqrt((p.coordinates(1)^2+p.coordinates(2)^2));
        for i = 1:length(Rdr(1,:))-1 
            if (Rdr(1,i) <= r) && (r < Rdr(1,i+1))
                x = i;       
            end
        end    
        if x == -1
            Escaped_bounds.coordinates = [Escaped_bounds.coordinates;p.coordinates];
	Escaped_bounds.weight = [Escaped_bounds.weight;p.weight];
Escaped_bounds.r = [Escaped_bounds.r;(-1)];
disp('Error -1');
        else
            if p.coordinates(3) == 0
                Rdr(2,x) = Rdr(2,x)+p.weight;
            elseif p.coordinates(3) == sum(tissue.layers(1:end))
                Tdr(2,x) = Tdr(2,x)+p.weight;
		else
		
            Escaped_bounds.coordinates = [Escaped_bounds.coordinates;p.coordinates];
	Escaped_bounds.weight = [Escaped_bounds.weight;p.weight];
Escaped_bounds.r = [Escaped_bounds.r;(2)];
disp('Error 2');
            end
        end        
        % 2D reflectance and transmittance arrays
        Rx = max_radius + 1; Ry = max_radius + 1; x = -1; y = -1;
        for i = 1:length(b) % for each available x coordinate
            if (p.coordinates(1)-b(i))^2 < Rx
                Rx = (p.coordinates(1)-b(i))^2;
                x = i; % index in coord system
            end
            if (p.coordinates(2)-b(i))^2 < Ry
                Ry = (p.coordinates(2)-b(i))^2;
                y = i; % index in coord system
            end
        end                
        if x == -1 || y == -1
disp('Error 2D');
            if p.layer_number == 0 % Reflectance
                R_e.c=[R_e.c;p.coordinates];
                R_e.p=[R_e.p;p.weight];
            elseif p.layer_number == length(tissue.layers)+1 % Transmittance
                T_e.c=[T_e.c;p.coordinates];
                T_e.p=[T_e.p;p.weight];
            else
                disp('Error 3');
            end
        else % Valid coordinates
            if p.layer_number == 0 % Reflectance
                if p.event_number == 0 % Unscattered       
                    R_unscat = R_unscat + p.weight;
                else % Diffuse reflectance
                    Rd(x,y,2) = Rd(x,y,2) + p.weight;
                end
            elseif p.layer_number == length(tissue.layers)+1 % Transmittance
                if p.event_number == 0 
                    T_unscat = T_unscat + p.weight;
                else
                    Td(x,y,2) = Td(x,y,2) + p.weight;
                end                
            else
                disp('Error 4');  
            end
        end
   end
end
% close(h); 
end
