function [] = fbm_print_table(results,scenario)
% call f.x. as fbm_print_table(load_super.results,1)
% where load_super=load('superdiffusive_data_and_results/superdiffusive_track_output.mat')

MM = [0 0 0 0; 1 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0; 1 1 0 1; 0 0 1 1; 1 1 1 1];
std_vals = {'$0$','$0$','$0$','$1/2$'};

if scenario==1
  misc.data_id = 'superdiffusive_data_and_results/superdiffusive_track';
  true_vals=[20 0 0 10 0.75];
  data_path = [misc.data_id,'.txt'];
  %loading data 
  data = load(data_path);
  data = transpose(data);
  %Subtract points to obtain steps
  obs = data(:,2:end) - data(:,1:end-1);
  trueLmax=fbm_logl(obs,true_vals,[1 1 1 1])/log(10);
  trueLstr=sprintf('%.2f',trueLmax)
elseif scenario==2
  misc.data_id = 'subdiffusive_data_and_results/subdiffusive_track';
  true_vals=[20 0 10 10 0.25];
  data_path = [misc.data_id,'.txt'];
  %loading data 
  data = load(data_path);
  data = transpose(data);
  %Subtract points to obtain steps
  obs = data(:,2:end) - data(:,1:end-1);
  trueLmax=fbm_logl(obs,true_vals,[1 1 1 1])/log(10);
  trueLstr=sprintf('%.2f',trueLmax)
else
  trueLstr=sprintf('-')
end

fid = fopen('print_table_output.txt','w');

for i=1:length(results)
  fprintf(fid,'%i & ',i);
  fprintf(fid,'$%.2f \\pm %.2f$ & ',results(i).logZ(1)/log(10),results(i).logZ_error/log(10));
  fprintf(fid,'$%.1f \\pm %.1f$ & ',results(i).param_mean(1),results(i).param_stddev(1));
  n=1;
  for j=1:length(MM(i,:))
    if MM(i,j)==1
      n=n+1;
      if j==4
        fprintf(fid,'$%.3f \\pm %.3f$ & ',results(i).param_mean(n),results(i).param_stddev(n));
      else
        fprintf(fid,'$%.1f \\pm %.1f$ & ',results(i).param_mean(n),results(i).param_stddev(n));
      end
    else
      fprintf(fid,[std_vals{j} ' & ']);
    end
  end
  fprintf(fid,'$%.2f$ & ',results(i).samples(end).logl/log(10));
  fprintf(fid,'$%.8f$ \\\\\n',results(i).Z_norm);
end

