function Q_EQ=isoterma_sips(PSI4,PSI3,n,X)
X(X<=0)=0;
Q_EQ=PSI4*PSI3*X.^(1/n)./(1+PSI3*X.^(1/n));
Q_EQ(Q_EQ<=0)=0;
end