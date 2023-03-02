% xls name
outputPubFolder = 'C:\Users\mgeva\Documents\GitHub\closedLoop-pub\figureGen';
sourceXls = fullfile(outputPubFolder,'sourceFig1.xls');

list = dir(fullfile(outputPubFolder,'*table*'));

%% Fig 1 - all panels with statistics included
sourceXls = fullfile(outputPubFolder,'sourceFig1.xls');

mm = matfile(fullfile(outputPubFolder,'Fig1c3_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','c3')

mm = matfile(fullfile(outputPubFolder,'Fig1f_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','f')

mm = matfile(fullfile(outputPubFolder,'Fig1g_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','g')


%% Fig 2 - all panels with statistics included
sourceXls = fullfile(outputPubFolder,'sourceFig2.xls');
mm = matfile(fullfile(outputPubFolder,'Fig2c_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','c')

mm = matfile(fullfile(outputPubFolder,'Fig2e1_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','e_1')

mm = matfile(fullfile(outputPubFolder,'Fig2e3_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','e_3')

mm = matfile(fullfile(outputPubFolder,'Fig2e2_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','e_2')

mm = matfile(fullfile(outputPubFolder,'Fig2f_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','f')

%% Fig 3 - all panels with statistics included
sourceXls = fullfile(outputPubFolder,'sourceFig3.xls');
mm = matfile(fullfile(outputPubFolder,'Fig3c_table.mat'));
T = mm.Fig3c_data_table;
writetable(T,sourceXls,'Sheet','c')

mm = matfile(fullfile(outputPubFolder,'Fig3d_table.mat'));
T = mm.Fig3d_data_table;
writetable(T,sourceXls,'Sheet','d')

mm = matfile(fullfile(outputPubFolder,'Fig3e_table.mat'));
T = mm.Fig3e_data_table;
writetable(T,sourceXls,'Sheet','e')

%% Fig 4 - all panels with statistics included
sourceXls = fullfile(outputPubFolder,'sourceFig4.xls');

mm = matfile(fullfile(outputPubFolder,'Fig4a_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','a')

mm = matfile(fullfile(outputPubFolder,'DataFig4c_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','c')

mm = matfile(fullfile(outputPubFolder,'DataFig4d_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','d')


%% Extended Data Fig 1 - 
sourceXls = fullfile(outputPubFolder,'sourceExtendedDataFig1.xls');

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig1c1_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','c1')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig1c2_table.mat'));
T = mm.jointTable;
writetable(T,sourceXls,'Sheet','c2')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig1c3_table.mat'));
T = mm.jointTable;
writetable(T,sourceXls,'Sheet','c3')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig1c4_table.mat'));
T = mm.jointTable;
writetable(T,sourceXls,'Sheet','c4')

%% Extended Data Fig 3 - 
sourceXls = fullfile(outputPubFolder,'sourceExtendedDataFig3.xls');

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig3c_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','c')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig3d_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','d')

%% Extended Data Fig 4 - 
sourceXls = fullfile(outputPubFolder,'sourceExtendedDataFig4.xls');

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig4a_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','a')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig4c_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','c')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig4d_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','d')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig4e_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','e')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig4g_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','g')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig4f_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','f')


%% Extended Data Fig 5 - 
sourceXls = fullfile(outputPubFolder,'sourceExtendedDataFig5.xls');

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig5a1_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','a')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig5a2_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','b')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig5a3_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','c')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig5_2_1_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','d')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig5_2_2_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','e')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig5_2_3_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','f')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig5_3_1_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','g')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig5_3_2_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','h')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig5_3_3_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','i')

%% Extended Data Fig 7 - 
sourceXls = fullfile(outputPubFolder,'sourceExtendedDataFig7.xls');
mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig1c1_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','a')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig7b_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','b')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig7c_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','c')



%% Extended Data Fig 9 - 
sourceXls = fullfile(outputPubFolder,'sourceExtendedDataFig9.xls');

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig9a1_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','a1')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig9a2_data_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','a2')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig9e_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','e')

mm = matfile(fullfile(outputPubFolder,'ExtendedDataFig9d_table.mat'));
T = mm.data_table;
writetable(T,sourceXls,'Sheet','d')

