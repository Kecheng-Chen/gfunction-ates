function extraction(alpha_t,alpha_Q,Q,lambda,rhoc,H,name_org)
    value=[alpha_t./(alpha_t+1),alpha_Q,Q/2,lambda,rhoc,H/2];
    [h1,h2]=extraction_all(value,name_org);
    figg3=figure(3);hold on;
    plot(1:50,h1);
    plot(1:50,h2);
    xlabel('Year');
    ylabel('h');
    legend('Injection process','Extraction process');
    title('Heat transfer coefficient during injection and extraction')
    saveas(figg3,strcat("results/heat_loss_coefficient.png"));
end