function fbm_skel_artificial()

%misc.data_id = 'superdiffusive_data_and_results/superdiffusive_track';
misc.data_id = 'subdiffusive_data_and_results/subdiffusive_track';

data_path = [misc.data_id,'.txt'];

%loading data 
data = load(data_path);
data = transpose(data);
%Subtract points to obtain steps
obs = data(:,2:end) - data(:,1:end-1);

logl_true = fbm_logl(obs,[20 0 0 10 0.75],[1 1 1 1])

%Tell ns_processdataset to write a summary file
misc.nssummary=['_results.txt'];

%Specify framerate and pixelsize
p_size = 1;  %Pixelsize in data in micrometer/pixel
tau = 0.5;     %Time between data point in seconds

%Specify prior ranges
sigmaHmin=10^(0);  %Minimum value for diffusion constant in micrometer^2/s
sigmaHmax=10^(3);      %Maximum value for diffusion constant in micrometer^2/s

vmax=1000;      %Maximum absolute value of drift in micrometer/s

noisemin=0;  %Minimum value for noise parameter in micrometers
noisemax=1000;  %Maximum value for noise parameter in micrometers

Hmin=0;      %Minimum value for Hurst parameter
Hmax=1;      %Maximum value for Hurst parameter

%Convert prior ranges to dimensions of pixels and frames for convenience 
% ranges - a 2x5 array of minimum and maximum values for the 5 parameter:
%   - diffusion coefficient, drift velocity (in two dimentions), noise
%     parameter, and Hurst exponent.
ranges=[sigmaHmin sigmaHmax;...
    -vmax vmax ;-vmax vmax;...
    noisemin noisemax; Hmin Hmax];

%Specify options
%options.nwalkers=10; % Number of walkers to be generated
%options.stoprat=10^(-3);
%options.nsteps=10;
options.nlist= [1 2 4 8 16 32];
options.trackmax = 100;

%Specify the models
final=@(x) x(end); % extract last element in a vector
MM = [0 0 0 0; 1 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0; 1 1 0 1; 0 0 1 1; 1 1 1 1];
for i=1:8;
  models(i).genu=@() util_generate_u(sum(MM(i,:))+1);
  models(i).options=options;
  models(i).logl=@(obs,theta) fbm_logl(obs,fbm_params(theta,MM(i,:)),MM(i,:));
  models(i).invprior=@(u) fbm_invprior(u,ranges,MM(i,:));
  models(i).scaling = @(obs,n) fbm_scaling(obs,n);
  models(i).replicate = @(obs,theta,n) fbm_replicate(obs,fbm_params(theta,MM(i,:)),n);
  models(i).logl_n = @(obs,theta,n) fbm_logl_n(obs,fbm_params(theta,MM(i,:)),MM(i,:),n);
  models(i).labels=[1];
  for j=1:length(MM(i,:));
    if MM(i,j)==1
      models(i).labels=[models(i).labels j+1];
    end
  end
  models(i).add{1}=@(theta) theta(1)^2/(2*tau^(2*final(fbm_params(theta,MM(i,:)))));
  models(i).labels=[models(i).labels 6];
  for j=1:2
    if MM(i,j)==1
      models(i).add{end+1}=@(theta) theta(j+1)/tau;
      models(i).labels=[models(i).labels 6+j];
    end
  end
end

%Labels for the parameters
misc.labels=...
['Step deviation: ';...
 'x-bias/step:    ';...
 'y-bias/step:    ';...
 'Measurem. err.: ';...
 'Hurst exponent: ';...
 'D_H constant:   ';...
 'x-velocity:     ';...
 'y-velocity:     '];

path=[misc.data_id,'_output'];

[results] = ns_processdataset(obs,models,misc);
save(path,'results')

