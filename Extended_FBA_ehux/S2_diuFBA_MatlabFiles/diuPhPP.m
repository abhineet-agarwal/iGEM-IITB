global diuFBA;
stepNo=50;

MultiPhenotypePhasePlane(1,diuFBA.model.cbmod,'t1_LIPID.ER','t1_mnl[c]_transfer',stepNo,.8,2);
MultiPhenotypePhasePlane(2,diuFBA.model.cbmod,'t1_protein.biomass1.6[c]_transfer','t1_mnl[c]_transfer',stepNo,-2,2);

     
