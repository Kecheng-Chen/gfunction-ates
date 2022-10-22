clear;clc;
% id of the example case study [1,2,3]
case_index=1;

% get parameters for case study
% alpha_t: ratio between the duration of the injection period and the
% duration of the extraction period. (dimensionless)
% alpha_Q: ratio between the injection flow rate and the extraction flow
% rate (dimensionless)
% Q: the base flow rate (m^3/s)
% lambda: heat conductivity (W/mK)
% rhoc: volumetric heat capacity (MJ/m^3K)
% H: aquifer thickness (m)
% name_org: the name of the file storing COMSOL generated benchmark results
[alpha_t,alpha_Q,Q,lambda,rhoc,H,name_org] = case_study(case_index);

% generate the comparison results
% saved in folder 'results/'
if alpha_t*alpha_Q>1*1
    % injection-dominated analytical scheme
    disp('Injection-dominated');
    injection(alpha_t,alpha_Q,Q,lambda,rhoc,H,name_org);
else
    % extraction-dominated analytical scheme
    disp('Extraction-doinated');
    extraction(alpha_t,alpha_Q,Q,lambda,rhoc,H,name_org);
end