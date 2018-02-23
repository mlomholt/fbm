Ntracks=170;
N=200; % generate 201 (x,y)-positions, i.e., 200 steps
tr_pos=zeros(N+1,2);
% Parameters in square brackets below: [sigma_H mu_x mu_y sigma_mn H]
MM = [0 0 0 0; 1 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0; 1 1 0 1; 0 0 1 1; 1 1 1 1];
model=4;

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

obslist={};
paramlist={};
modellist={};

for i=1:Ntracks
  model=randi(8);
  u=util_generate_u(sum(MM(model,:))+1);
  theta=fbm_invprior(u,ranges,MM(model,:));
  params=fbm_params(theta,MM(model,:));
  tr_steps=fbm_replicate(zeros(N,2),params,1); % Generate steps
  tr_pos(2:(N+1),:)=(cumsum(transpose(tr_steps))); % Calculate positions
  obslist{i}=tr_pos;
  paramlist{i}=params;
  modellist{i}=model;
end

save('random_data_and_results/random_tracks.mat','obslist','paramlist','modellist')
