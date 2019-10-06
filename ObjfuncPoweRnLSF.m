function [ff] = ObjfuncPoweRnLSF(TPL,LSF)

W1 = 0.8;
W2 =0.2;
% 
% DistLoadFlowSolution=powerflow;
% 
% VTM=DistLoadFlowSolution.VmagPU;
% 
% TPL=DistLoadFlowSolution.PtLosskW;

ff = W1*TPL + W2*sum((1-LSF).^2);         %TPL = TotalPowerLoss;  VTM=VoltageMag
% ff = 1./(1+(f.^2));
end

