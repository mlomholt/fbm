function [] = fbm_print_confusionmatrix(confusionmat)
% call f.x. as fbm_print_table(load_super.results,1)
% where load_super=load('superdiffusive_data_and_results/superdiffusive_track_output.mat')

fid = fopen('print_confusionmatrix.txt','w');

fprintf(fid,'%% Row numbering represent true models\n');
fprintf(fid,'%% Column numbering represent infered models\n');
fprintf(fid,'  ');
for j=1:length(confusionmat(1,:))
  fprintf(fid,'& %i ',j);
end
for i=1:length(confusionmat(:,1))
fprintf(fid,'\\\\\n');
  fprintf(fid,'%i ',i);
  for j=1:length(confusionmat(1,:))
    fprintf(fid,'& %i ',confusionmat(i,j));
  end
end
fprintf(fid,'\n');
fclose(fid);

