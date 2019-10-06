 
%% This file will print all Load-flow results

% Load-flow results are also stored in DistLoadFlowSolution
clear all
clc
DistLoadFlowSolution=powerflow
display('|| All Results are stored in DistLoadFlowSolution (see in Work-space)||')
display('|| Some Results, such as, Voltage and Losses are also displayed here ||')
Vsol=[[1:length(DistLoadFlowSolution.Vactual)]' DistLoadFlowSolution.Vactual];
display('Voltage Solution in kV is displayed below...')
display('    Bus No.  Voltage(kV)')
disp(Vsol)
Losses=[[1:length(DistLoadFlowSolution.Pbrloss)]' DistLoadFlowSolution.Pbrloss DistLoadFlowSolution.Qbrloss];
display('Losses in each Line segment are displayed in kWs and kVAr...')
display('   Line No.  Loss(kWs) Loss(kVAr)')
disp(Losses)
disp(' Total Lines Losses in kW, kVAr, and in kVA are...')
display('  TotalkW   TotalkVAr  TotalkVA')
disp([DistLoadFlowSolution.PtLosskW DistLoadFlowSolution.QtLosskVAr DistLoadFlowSolution.SLosskVA])
clearvars -except DistLoadFlowSolution