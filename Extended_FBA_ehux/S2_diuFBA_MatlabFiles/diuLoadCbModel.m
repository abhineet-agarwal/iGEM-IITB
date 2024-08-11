function diuLoadCbModel()
% Load a Cobra Toolbox model for diuFBA

global diuFBA;

%% Load Cobra Model

load iEH410 cbmod
diuFBA.integrals=[findMetIDs(cbmod,'mnl[c]') findMetIDs(cbmod,'lipid.biomass5.[c]') findMetIDs(cbmod,'protein.biomass1.6[c]') findMetIDs(cbmod,'lipid.chloroplast2[c]') findMetIDs(cbmod,'biomass_log_200_24h_neutrcharge_final[c]')]
diuFBA.objective=[zeros(length(cbmod.c),1);0;0;0;0;-1];
%cbmod.ub(190)=0; findMetIDs(cbmod,'13glucan(300)[c]') 
%cbmod.ub(191)=0;

%% Assign COBRA model elements to the according diuFBA.model elements

diuFBA.S=cbmod.S;
diuFBA.cbmod=cbmod;
diuFBA.mNames=cbmod.metNames;
diuFBA.rNames=cbmod.rxnNames;
diuFBA.mIDs=cbmod.mets;
diuFBA.rIDs=cbmod.rxns;
for k=1:length(diuFBA.integrals)
    diuFBA.rNames{k+size(diuFBA.S,2)}=[diuFBA.mNames{diuFBA.integrals(k)},'_transfer'];
    diuFBA.rIDs{k+size(diuFBA.S,2)}=[diuFBA.mIDs{diuFBA.integrals(k)},'_transfer'];
end

diuFBA.ubSingle=[cbmod.ub;1000*ones(length(diuFBA.integrals),1)];
diuFBA.lbSingle=[cbmod.lb;-1000*ones(length(diuFBA.integrals),1)];
diuFBA.rev=[cbmod.rev;zeros(length(diuFBA.integrals),1)];
diuFBA.rNum=length(diuFBA.rev);
diuFBA.mNum=length(diuFBA.mNames);

