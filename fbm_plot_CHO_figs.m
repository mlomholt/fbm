function [] = fbm_plot_CHO_figs(load_25,load_37,load_random,load_track)

%if scenario==-1
%  misc.data_id = 'superdiffusive_data_and_results/superdiffusive_track';
%  data_path = [misc.data_id,'.txt'];
%  %loading data 
%  data = load(data_path);
%  data = transpose(data);
%  %Subtract points to obtain steps
%  obs = data(:,2:end) - data(:,1:end-1);
%elseif scenario==-2
%  misc.data_id = 'subdiffusive_data_and_results/subdiffusive_track';
%  data_path = [misc.data_id,'.txt'];
%  %loading data 
%  data = load(data_path);
%  data = transpose(data);
%  %Subtract points to obtain steps
%  obs = data(:,2:end) - data(:,1:end-1);
%elseif scenario==1

%load_25=load('CHO_data_and_results/CHO_BC_25_2_output.mat')
[countsi1,imaxlist1,Hlist1,Hlisttrue1,Dlist1,Dlisttrue1,pinflist1] = analyze(load_25.res_frederik);
%load_37=load('CHO_data_and_results/CHO_BC_37_output.mat')
[countsi2,imaxlist2,Hlist2,Hlisttrue2,Dlist2,Dlisttrue2,pinflist2] = analyze(load_37.res_frederik);
%load_random=load('random_data_and_results/random_tracks_output.mat');
[countsi3,imaxlist3,Hlist3,Hlisttrue3,Dlist3,Dlisttrue3,pinflist3] = analyze(load_random.res_frederik);

misc.data_id = 'CHO_data_and_results/CHO_BC_25_2';
data_path = [misc.data_id,'.mat'];
data = load(data_path);  %loading data 
%Extract track coordinates
Trajectory1=1;
xtraj1=data.x25_2(:,Trajectory1);
ytraj1=data.y25_2(:,Trajectory1);
B_2d1 = transpose([xtraj1 ytraj1]);

H4vals=zeros(1,length(load_25.res_frederik));
for i=1:length(H4vals)
  H4vals(i)=load_25.res_frederik{i}(4).maxLpar(2);
end
[H,Trajectory2]=max(H4vals)
mod4Q=imaxlist1(Trajectory2)
xtraj2=data.x25_2(:,Trajectory2);
ytraj2=data.y25_2(:,Trajectory2);
B_2d2 = transpose([xtraj2 ytraj2]);
 

countdist=zeros(length(load_random.res_frederik{1}),2);
truedist=zeros(length(load_random.res_frederik{1}),2);
for i=1:length(load_track.modellist)
  if load_track.modellist{i}==imaxlist3(i)
   countdist(imaxlist3(i),1)=countdist(imaxlist3(i),1)+1;
   truedist(load_track.modellist{i},1)=truedist(load_track.modellist{i},1)+1;
  else
   countdist(imaxlist3(i),2)=countdist(imaxlist3(i),2)+1;
   truedist(load_track.modellist{i},2)=truedist(load_track.modellist{i},2)+1;
  end
end

correct=sum(countdist(:,1))
wrong=sum(countdist(:,2))

figure(1)
subplot(2,2,1)
plot(xtraj1,ytraj1,'-')
subplot(2,2,2)
plot(xtraj2,ytraj2,'-')
subplot(2,2,3)
histogram(Hlist1,'BinLimits',[0 1])
subplot(2,2,4)
histogram(Hlist2,'BinLimits',[0 1])
%bar(1:8,countsi3)

figure(2)
subplot(2,2,1)
bar(1:8,countsi1)
subplot(2,2,2)
bar(1:8,countsi2)
subplot(2,2,3)
bar(countdist,'stacked')
%histogram(Hlist3,'BinLimits',[0 1])
subplot(2,2,4)
bar(truedist,'stacked')
%histogram(Hlisttrue3)
%histogram(imaxlist,'BinEdges',0.5:1:8.5)

figure(3)
subplot(2,2,1)
histogram(pinflist3{1})
subplot(2,2,2)
histogram(pinflist3{2})
subplot(2,2,3)
histogram(pinflist3{3})
subplot(2,2,4)
histogram(pinflist3{5})

figure(4)
subplot(2,2,1)
histogram(pinflist1{1},'BinLimits',[0 1])
subplot(2,2,2)
histogram(pinflist1{2},'BinLimits',[0 1])
subplot(2,2,3)
histogram(pinflist1{3},'BinLimits',[0 1])
subplot(2,2,4)
histogram(pinflist1{5},'BinLimits',[0 1])

%---------------------------------
function [countsi,imaxlist,Hlist,Hlisttrue,Dlist,Dlisttrue,pinflist] = analyze(results)

imaxlist=[];

modH=4;
paramH=2;
Hlist=[];
Hlisttrue=[];

modD=1;
paramD=1;
Dlist=[];
Dlisttrue=[];

for i=1:length(results)
  logZlist=[];
  for j=1:length(results{1})
    logZlist=[logZlist results{i}(j).logZ(1)];
  end
  [Zmax imax]=max(logZlist);
  imaxlist=[imaxlist imax];
  Hlist=[Hlist results{i}(modH).maxLpar(paramH)];
  if imax==modH %| imax==1
    Hlisttrue=[Hlisttrue results{i}(modH).maxLpar(paramH)];
  end
  Dlist=[Dlist results{i}(modD).maxLpar(paramD)];
  if imax==modD
    Dlisttrue=[Dlisttrue results{i}(modD).maxLpar(paramD)];
  end
end
countsi=zeros(1,length(results{1}));
for i=1:length(countsi)
  countsi(i)=sum(imaxlist==i);
end

pinflist={};
for probi=1:length(results{i}(1).prob)
  pinflist{probi}=[];
  for i=1:length(results)
    pinflist{probi}=[pinflist{probi} results{i}(imaxlist(i)).prob(probi)];
  end
end

