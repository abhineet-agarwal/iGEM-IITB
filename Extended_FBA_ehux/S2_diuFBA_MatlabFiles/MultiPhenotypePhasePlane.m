function [growthRates,shadowPrices1,shadowPrices2] = MultiPhenotypePhasePlane(PlotNr,model,controlRxn1,controlRxn2,nPts,range1,range2,ShadowPlots,ShadowMet1,ShadowMet2)
%phenotypePhasePlane Plots three phenotype phase planes for two reactions.  The first plot is 
% a double robustness analysis, a kind of 3D surface plot.  The second
% two plots show the shadow prices of the metabolites from the two control 
% reactions, which define the phases.  Use the COLORMAP and SHADING 
% functions to change the looks of the plots.
%
% [growthRates,shadowPrices1,shadowPrices2] = phenotypePhasePlane(model,controlRxn1,controlRxn2,nPts,range1,range2)
%
%INPUTS
% PlotNr            Number of plot(s) in a series of plots (DK20150302)
% model             COBRA model structure
% controlRxn1       the first reaction to be plotted
% controlRxn2       the second reaction to be plotted
%
%OPTIONAL INPUTS
% nPts              the number of points to plot in each dimension
%                   (Default = 50)
% range1            the range of reaction 1 to plot
%                   (Default = 20)
% range2            the range of reaction 2 to plot
%                   (Default = 20)
% ShadowPlots       Plot the shadow prices for reaction 1 and 2?
%                   (Default = false) (DK20150302)
%
%OUTPUTS
% growthRates1      a matrix of maximum growth rates
% shadowPrices1     a matrix of rxn 1 shadow prices
% shadowPrices2     a matrix of rxn 2 shadow prices
%
% Jeff Orth 6/26/08
% David Knies 3/2/15

if (nargin < 5)
    nPts = 50;
end
if (nargin < 6)
    range1 = 20;
end
if (nargin < 7)
    range2 = 20;
end
if (nargin < 8)
    ShadowPlots=false;
end


% find rxn and met ID numbers to get shadow prices and reduced costs
rxnID1 = findRxnIDs(model,controlRxn1);
metID1 = find(model.S(:,rxnID1));
rxnID2 = findRxnIDs(model,controlRxn2);
metID2 = find(model.S(:,rxnID2));

if (nargin < 9)
    if ShadowPlots==true
        disp('Please select metabolite 1 for Shadow Price Plot:')
        for k=1:length(metID1)
            fprintf('%i\t%s\n',k,model.mets{metID1(k)});
        end
        ShadowMet1=input('Metabolite 1: ');
    end
end
if (nargin < 10)
    if ShadowPlots==true
        disp('Please select metabolite 2 for Shadow Price Plot:')
        for k=1:length(metID2)
            fprintf('%i\t%s\n',k,model.mets{metID2(k)});
        end
        ShadowMet2=input('Metabolite 2: ');
    end
end

        


% create empty vectors for the results
ind1 = linspace(0,range1,nPts);
ind2 = linspace(0,range2,nPts);
growthRates = zeros(nPts);
shadowPrices1 = zeros(nPts);
shadowPrices2 = zeros(nPts);

% calulate points
h = waitbar(0,sprintf('generating PhPP No. %i',PlotNr));
global CBT_LP_PARAMS  % save the state of the primal only flag.  
if isfield( CBT_LP_PARAMS, 'primalOnly')
    primalOnlySave = CBT_LP_PARAMS.primalOnly;
end
changeCobraSolverParams('LP', 'primalOnly', false);
for i = 1:nPts %ind1
    for j = 1:nPts %ind2
        waitbar((nPts*(i-1)+j)/(nPts^2),h);
        % (DK20150302) Changed sign to allow arbitrary reactions
        model1 = changeRxnBounds(model,controlRxn1,1*ind1(i),'b');
        model1 = changeRxnBounds(model1,controlRxn2,1*ind2(j),'b');
                
        fbasol = optimizeCbModel(model1,'max',0,true);
        if isfield(fbasol,'vbasis')
            model.cbasis=fbasol.cbasis;
            model.vbasis=fbasol.vbasis;
        end
        try
            growthRates(j,i)= fbasol.f;
            % to plot other fluxes on the z axis:
            % growthRates(j,i) = max(abs([fbasol.x(findRxnIDs(model,'t1_EX_co2(e)')),fbasol.x(findRxnIDs(model,'t2_EX_co2(e)'))]));
        end
        try % calculate shadow prices
            shadowPrices1(j,i) = fbasol.y(metID1(ShadowMet1));
            shadowPrices2(j,i) = fbasol.y(metID2(ShadowMet2));
        end
    end
end
if ( regexp( version, 'R20') )
        close(h);
end
if exist('primalOnlySave', 'var')
    changeCobraSolverParams('LP', 'primalOnly', primalOnlySave);
else
    changeCobraSolverParams('LP', 'primalOnly', true);
end

% % plot the points

% (DK20150302) Plot shadow prices? Adapted labels for diuFBA
if(ShadowPlots)
figure((PlotNr-1)*3+2);
pcolor(ind1,ind2,shadowPrices1);
xlabel(strrep(strcat(controlRxn1,' (mmol/g DW-hr)'),'_','\_')), ylabel(strrep(strcat(controlRxn2,' (mmol/g DW-hr)'),'_','\_')), zlabel('growth rate (1/d)');
figure((PlotNr-1)*3+3);
pcolor(ind1,ind2,shadowPrices2);
xlabel(strrep(strcat(controlRxn1,' (mmol/g DW-hr)'),'_','\_')), ylabel(strrep(strcat(controlRxn2,' (mmol/g DW-hr)'),'_','\_')), zlabel('growth rate (1/d)');
figure((PlotNr-1)*3+1);
surfl(ind1,ind2,growthRates);
xlabel(strrep(strcat(controlRxn1,' (mmol/g DW-hr)'),'_','\_')), ylabel(strrep(strcat(controlRxn2,' (mmol/g DW-hr)'),'_','\_')), zlabel('growth rate (1/d)');
figure(PlotNr+20)
surf(ind1,ind2,growthRates);
else
% (DK20150302) added PhPP only, adapted labels for diuFBA    
figure((PlotNr-1)*2+1);
% Plot with plane angle based shading
surfl(ind1,ind2,growthRates);
xlabel(strrep(strcat(controlRxn1,' (mmol/g DW)'),'_','\_')), ylabel(strrep(strcat(controlRxn2,' (mmol/g DW)'),'_','\_')), zlabel('growth rate (1/d)');
figure((PlotNr-1)*2+2)
% Plot with colormap
surf(ind1,ind2,growthRates);
xlabel(strrep(strcat(controlRxn1,' (mmol/g DW)'),'_','\_')), ylabel(strrep(strcat(controlRxn2,' (mmol/g DW)'),'_','\_')), zlabel('growth rate (1/d)');
end


