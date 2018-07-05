Ntracks=111;
%N=199; % generate N+1 (x,y)-positions, i.e., N steps
N=99; % generate N+1 (x,y)-positions, i.e., N steps
tr_pos=zeros(N+1,2);
% Parameters in square brackets below: [sigma_H mu_x mu_y sigma_mn H] (if applicable
  model=4;
  MMp=[0 0 0 1];

obslist={};
paramlist={};
modellist={};

for i=1:Ntracks
  theta=[10 0.6];
  params=fbm_params(theta,MMp);
  tr_steps=fbm_replicate(zeros(N,2),params,1); % Generate steps
  tr_pos(2:(N+1),:)=(cumsum(transpose(tr_steps))); % Calculate positions
  obslist{i}=tr_pos;
  paramlist{i}=params;
  modellist{i}=model;
end

%save('random_data_and_results/equivalent_tracks.mat','obslist','paramlist','modellist')
save('random_data_and_results/equivalent_short_tracks.mat','obslist','paramlist','modellist')
