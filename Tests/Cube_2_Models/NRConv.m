

Elastic_rhs=[1.000000e+00 7.983962e-03 2.763450e-05 3.243472e-08 5.078036e-11];
Elastic_iter=0:length(Elastic_rhs)-1;

Plastic_rhs=[1.000000e+00 4.245456e-01 1.835625e-01 2.313564e-01 2.928492e-01 7.046108e-02 2.588838e-04 5.480182e-07 2.496837e-09];
Plastic_iter=0:length(Plastic_rhs)-1;

figure (9)
semilogy(Elastic_iter,Elastic_rhs,'b-*',LineWidth=1.5)
hold on
semilogy(Plastic_iter,Plastic_rhs,'r-*',LineWidth=1.5)
xlabel('Iterazioni','FontSize',15)
ylabel('|rhs|/|rhs_0|','FontSize',15)
title('Profilo di convergenza NR globale','FontSize',18)
legend('Fase elastica','Fase plastica','FontSize',15)
grid on
