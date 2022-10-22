function [alpha_t,alpha_Q,Q,lambda,rhoc,H,name_org] = case_study(case_index)
switch case_index
    case 1
        alpha_t=0.55./0.45;
        alpha_Q=1.2;
        Q=0.02;
        lambda=2.8;
        rhoc=2;
        H=16;
        name_org='validation_data/test-sol1.txt';
    case 2
        alpha_t=0.6./0.4;
        alpha_Q=2;
        Q=0.0417;
        lambda=3.01;
        rhoc=2.74;
        H=40;
        name_org='validation_data/test-sol2.txt';
    case 3
        alpha_t=0.5./0.5;
        alpha_Q=0.5;
        Q=0.0056;
        lambda=2.51;
        rhoc=1.64;
        H=10;
        name_org='validation_data/test-sol3.txt';
end
end