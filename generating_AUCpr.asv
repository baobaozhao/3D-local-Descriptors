function Area = generating_AUCpr(P,R)
% generating AUCpr value based on PRC.

P=1-P;
L=length(R);
Area=0;
for i=1:L-1
    if i==1
    X_zhi=R(i+1)-0;
    else
    X_zhi=R(i+1)-R(i); 
    end
    Y_zhi=(P(i+1)+P(i))/2;
    Area=Area+X_zhi*Y_zhi;
    
end


end

