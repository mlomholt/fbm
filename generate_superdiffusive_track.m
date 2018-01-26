N=200; % generate 201 (x,y)-positions, i.e., 200 steps
tr_pos=zeros(N+1,2);
% Parameters in square brackets below: [sigma_H mu_x mu_y sigma_mn H]
tr_steps=fbm_replicate(zeros(N,2),[20 0 0 10 3/4],1); % Generate steps
tr_pos(2:(N+1),:)=(cumsum(transpose(tr_steps))); % Calculate positions
save('superdiffusive_data_and_results/superdiffusive_track.txt','tr_pos','-ascii','-tabs')
