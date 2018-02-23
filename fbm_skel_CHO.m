function fbm_skel_CHO()

%misc.data_id = 'CHO_data_and_results/CHO_BC_25_2'; Ntraj=111;
%misc.data_id = 'CHO_data_and_results/CHO_BC_37'; Ntraj=170;
misc.data_id = 'random_data_and_results/random_tracks'; Ntraj=170;

data_path = [misc.data_id,'.mat'];

%loading data 
data = load(data_path);
obs={};
for i=1:Ntraj
  %Extract track coordinates
  if misc.data_id(end)=='2'
    xtraj=1000*data.x25_2(:,i);
    ytraj=1000*data.y25_2(:,i);
  elseif misc.data_id(end)=='7'
    xtraj=1000*data.x(:,i);
    ytraj=1000*data.y(:,i);
  else
    xtraj=data.obslist{i}(:,1);
    ytraj=data.obslist{i}(:,2);
  end
  B_2d = transpose([xtraj ytraj]);

  %transpose([x y]);
  %Subtract points to obtain steps
  X_2d= B_2d(:,2:end) - B_2d(:,1:end-1);
  obs=horzcat(obs,X_2d);
end

lenobs=length(obs);

%Tell ns_processdataset to write a summary file
misc.nssummary=['_results.txt'];

%Specify framerate and pixelsize
p_size = 1000;  %Pixelsize in data in nanometer/pixel
tau = 0.5;     %Time between data point in seconds

%Specify prior ranges
sigmaHmin=10^(0);  %Minimum valu for step deviation in nanometers
sigmaHmax=10^(3);      %Maximum valu for step deviation in nanometers

vmax=1000;      %Maximum absolute value of drift in nanometer/s

noisemin=0;  %Minimum value for noise parameter in nanometers
noisemax=1000;  %Maximum value for noise parameter in nanometers

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
options.nlist= [1 2 4 8 16 32];
options.trackmax = 100;
options.maxsamples=1000;

%Specify the models
final=@(x) x(end); % extract last element in a vector
MM = [0 0 0 0; 1 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0; 1 1 0 1; 0 0 1 1; 1 1 1 1];
for i=1:8;
  nparams=sum(MM(i,:))+1;
  models(i).genu=@() util_generate_u(nparams);
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

res_frederik={};
for i=1:lenobs
  misc.append = ['--- Trajectory ' num2str(i) ' ---\n\n'];
  [results] = ns_processdataset(obs{i},models,misc);
  res_frederik{i}=results;
  save(path,'res_frederik')
  fprintf('Individual track number %i analyzed\n',i);
end

