%% dFBA for not growth optimal Organisms
global diuFBA;

try
    if isempty(whos('global','CBT_LP_PARAMS'))
        initCobraToolbox
    end
catch
    disp('Cobra toolbox required')
end
%% Load external data and set parameters

diuLoadCbModel;
%diuFBALoadCNAModel;
diuFBA.timesteps=2;
%% Create extended stoichiometric matrix and reaction labels
tic;
[diuFBA.model.A,diuFBA.rLabel,diuFBA.mLabel,diuFBA.rIDs,diuFBA.mIDs]=diuSbuilder;
buildTime=toc;
disp(['Time for system matrix construction: ',num2str(buildTime),'s'])

diuFBA.crLabel=char(diuFBA.rLabel);
diuFBA.cmLabel=char(diuFBA.mLabel);
diuFBA.crIDs=char(diuFBA.rIDs);
diuFBA.cmIDs=char(diuFBA.mIDs);
%% Set up gurobi structure
% Create bounds for extended system
diuSetConstraints;
% Minimization coefficients: -1 for growth
diuFBA.model.sense='=';
% Objective function evaluated in last time step
diuFBA.model.obj=[zeros((diuFBA.timesteps-1)*diuFBA.rNum,1);diuFBA.objective];
diuFBA.model.rhs=zeros(size(diuFBA.model.A,1),1);

% Cobra model
diuFBA.model.cbmod.S=diuFBA.model.A;
diuFBA.model.cbmod.lb=diuFBA.model.lb;
diuFBA.model.cbmod.ub=diuFBA.model.ub;
diuFBA.model.cbmod.c=-diuFBA.model.obj;
diuFBA.model.cbmod.rxns=diuFBA.rIDs;
diuFBA.model.cbmod.mets=diuFBA.mIDs;
diuFBA.model.cbmod.rxnNames=diuFBA.rLabel;
diuFBA.model.cbmod.metNames=diuFBA.mLabel;
diuFBA.model.cbmod.rev=zeros(size(diuFBA.model.cbmod.lb));
diuSetConstraints;
for i=1:length(diuFBA.model.cbmod.lb)
    if diuFBA.model.cbmod.lb<0
        diuFBA.model.cbmod.rev(i)=1;
    end
end

%% direct Gurobi optimization
diuFBA.res.unref=gurobi(diuFBA.model);

%% Cobra Toolbox optimization
%diuFBA.res.unref=optimizeCbModel(diuFBA.model.cbmod);



fOutput=fopen('diuFBAres.txt','w');
for i=1:length(diuFBA.res.unref.x)
    disp([diuFBA.crIDs(i,:),' | ',num2str(diuFBA.res.unref.x(i),'%3.3f\n')])
    fprintf(fOutput,[diuFBA.crLabel(i,:),' | ',num2str(diuFBA.res.unref.x(i),'%3.3f'),'\n']);
    if(mod(i,diuFBA.rNum)==0)
        disp(' ');
        fprintf(' \n');
    end
end
fclose(fOutput);
for l=1:diuFBA.timesteps
    for k=1:diuFBA.rNum
       diuFBA.res.mat(k,l)=diuFBA.res.unref.x((l-1)*diuFBA.rNum+k);
    end
end
