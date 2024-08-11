function [Stilde,reac_label,met_label,reac_ID,met_ID]=diuSbuilder
global diuFBA;
tic;
% Create Transfer Matrix T
T=zeros(2*diuFBA.mNum,length(diuFBA.integrals));
for k=1:length(diuFBA.integrals)
    T(diuFBA.integrals(k),k)=-1;
    T(diuFBA.integrals(k)+diuFBA.mNum,k)=1;    
end
% Create block matrix S|T
blck=sparse([[diuFBA.S;zeros(size(diuFBA.S))],T]);
% initialize full S matrix Stilde
Stilde=blck;
% place blocks on "diagonal"
for k=1:diuFBA.timesteps-2
%    for l=1:(2*diuFBA.mNum)
%        for m=1:diuFBA.rNum
%            Sdyn((k-1)*diuFBA.mNum+l,(k-1)*diuFBA.rNum+m)=blck(l,m);
%        end
%    end  
    Stilde=[[Stilde;sparse(zeros(diuFBA.mNum,k*diuFBA.rNum))],[sparse(zeros((2*(k-1)+1)*diuFBA.mNum,diuFBA.rNum));blck]];
end
% last Block with incomplete T matrix
for l=1:(diuFBA.mNum)
   for m=1:size(diuFBA.S,2)
       Stilde((diuFBA.timesteps-1)*diuFBA.mNum+l,(diuFBA.timesteps-1)*diuFBA.rNum+m)=diuFBA.S(l,m);
   end
end
for l=1:(diuFBA.mNum)
   for m=1:length(diuFBA.integrals)
       Stilde((diuFBA.timesteps-1)*diuFBA.mNum+l,(diuFBA.timesteps-1)*diuFBA.rNum+m+size(diuFBA.S,2))=T(l,m);
   end
end

time1=toc;
% create reaction and metabolite label placeholders, to be replaced in main file
reac_label={};
reac_ID={};
met_label={};
met_ID={};
for k=1:diuFBA.timesteps
    for m=1:diuFBA.rNum
        reac_label{(k-1)*diuFBA.rNum+m}=['t',num2str(k),'_',diuFBA.rNames{m}];
        reac_ID{(k-1)*diuFBA.rNum+m}=['t',num2str(k),'_',diuFBA.rIDs{m}];
    end
end
for k=1:diuFBA.timesteps
    for m=1:diuFBA.mNum
        met_label{(k-1)*diuFBA.mNum+m}=['t',num2str(k),'_',diuFBA.mNames{m}];
        met_ID{(k-1)*diuFBA.mNum+m}=['t',num2str(k),'_',diuFBA.mIDs{m}];
    end
end

time2=toc;
fprintf('\n Time for matrix construction: %4.2f\n Time for label construction: %4.2f\n',time1,time2);

