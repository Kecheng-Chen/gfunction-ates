function injection(alpha_t,alpha_Q,Q,lambda,rhoc,H,name_org)
    value=[alpha_t,alpha_Q,Q/2,lambda,rhoc,H/2];
    bool=pyrunfile("sol.py","c",X=[value(1),value(2),value(3),value(4),value(5),value(6),1]);
    if bool==0
        [h1,h2]=injection_advection(value,name_org);
    else
        [h1,h2]=injection_diffusion(value,name_org);
    end
    figg3=figure(3);hold on;
    plot(1:50,h1);
    plot(1:50,h2);
    xlabel('Year');
    ylabel('h');
    legend('Injection process','Extraction process');
    title('Heat transfer coefficient during injection and extraction')
    saveas(figg3,strcat("results/heat_loss_coefficient.png"));
end