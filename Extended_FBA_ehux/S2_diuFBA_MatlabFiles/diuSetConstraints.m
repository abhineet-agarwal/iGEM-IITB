function [ output_args ] = diuSetConstraints(  )
%DIUSETCONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

global diuFBA;
diuFBA.model.lb=[];
diuFBA.model.ub=[];
for k=1:2            
    diuFBA.model.lb=[diuFBA.model.lb;diuFBA.lbSingle];
    diuFBA.model.ub=[diuFBA.model.ub;diuFBA.ubSingle];
end
    diuFBA.model.lb=diuFBA.model.lb*10;
    diuFBA.model.ub=diuFBA.model.ub*10;

for k=0:(length(diuFBA.integrals)-1)
    diuFBA.model.lb(length(diuFBA.model.lb)-k)=0;
end
   
if(isfield(diuFBA.model,'cbmod'))
diuFBA.model.ub(findRxnIDs(diuFBA.model.cbmod,'t1_GROWTH_log_200_16h_neutralcharge'))=0;
diuFBA.model.ub(findRxnIDs(diuFBA.model.cbmod,'t1_GROWTH_log_200_24h_neutralcharge'))=0;
diuFBA.model.ub(findRxnIDs(diuFBA.model.cbmod,'t1_GROWTH_log_50_16h_neutralcharge'))=0;
diuFBA.model.ub(findRxnIDs(diuFBA.model.cbmod,'t1_GROWTH_log_50_24h_neutralcharge'))=0;
diuFBA.model.lb(findRxnIDs(diuFBA.model.cbmod,'t1_CALCI-2-O'))=1.47*16;
diuFBA.model.lb(findRxnIDs(diuFBA.model.cbmod,'t2_CALCI-2-O'))=1.47*8;
diuFBA.model.lb(findRxnIDs(diuFBA.model.cbmod,'t1_EX_hn(e)'))=-24.24*16;
diuFBA.model.lb(findRxnIDs(diuFBA.model.cbmod,'t2_EX_hn(e)'))=0;

diuFBA.model.cbmod.lb=diuFBA.model.lb;
diuFBA.model.cbmod.ub=diuFBA.model.ub;

end


