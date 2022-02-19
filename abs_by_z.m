
% --- Extract absorption data by depth ---

function [] = abs_by_z(project_name)

% Find experiment numbers in each project
E = MCS_E(project_name); 

parameter_file = strcat('Data/Exp',sprintf('%05d',E(1)),'/parameters.mat');
if ~~exist(parameter_file)
load(parameter_file); 
z_bins = linspace(0,sum(parameters.tissue.layers(:)),101);
Z = z_bins(1:end-1);

A = zeros(length(E),100); nodes = zeros(2,length(E)); P = zeros(1,length(E));
for e = 1:length(E)
    filename = strcat('Data/Exp',sprintf('%05d',E(e)),'Exp',sprintf('%05d',E(e))); 
    z = []; a = []; photons = 0;
    for n = 1:parameters.simulation.number_of_nodes
        disp(strcat('___',num2str(n),'/',num2str(parameters.simulation.number_of_nodes)));
        f = strcat('Data/Exp',sprintf('%05d',E(e)),'/',num2str(n),'.mat');
        try load(f,'abs_coords');
            R = [abs_coords{:}]; clear('abs_coords');
            R = reshape(R,3,length(R)/3);
            z = [z;R(3,:)];
            photons = photons + parameters.simulation.number_of_photons;
        catch
            % Skip
            nodes(1,e) = 1;
        end
        try load(f,'abs_weight');
            aa = [abs_weight{:}]; clear('abs_weight');
            a = [a;aa];  
        catch
            % Skip
            nodes(2,e) = 1;
        end
    end
idx = zeros(length(z_bins)-1,length(z)); Abs = zeros(1,length(z_bins)-1);
if ~isempty(a)
for i = 1:length(z_bins)-1
    idx(i,:) = (z>=z_bins(i))&(z<z_bins(i+1)); 
    Abs(i) = sum(a(logical((idx(i,:)))),2);
end
end

A(e,:) = Abs;
P(e) = photons;
disp(strcat(num2str(e),'/',num2str(length(E))));

end
save(strcat('Results/',project_name{1},'_A_Z'),'A','Z','P','nodes');
end
end
