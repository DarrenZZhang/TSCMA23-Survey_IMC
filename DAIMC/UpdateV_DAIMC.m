function V = UpdateV_DAIMC(X,W,U,V,viewNum)

time = 0;
f = 0;
while 1
    time = time +1;
    sumVUUminus = 0;
    sumVUUplus = 0;
    sumXUminus = 0;
    sumXUplus = 0;
    for i = 1:viewNum
        XU = X{i}'*U{i};
        absXU = abs(XU);
        XUplus = (absXU + XU)/2;
        XUminus = (absXU - XU)/2;
        
        UU = U{i}'*U{i};
        absUU = abs(UU);
        UUplus = (absUU + UU)/2;
        UUminus = (absUU - UU)/2;
            
        sumXUminus = sumXUminus + W{i}*XUminus;
        sumXUplus = sumXUplus + W{i}*XUplus;
        
        sumVUUplus = sumVUUplus + W{i}*V*UUplus;
        sumVUUminus = sumVUUminus + W{i}*V*UUminus;
    end

    V = V.*sqrt((sumXUplus + sumVUUminus)./(max(sumXUminus + sumVUUplus,1e-10)));  
    
    ff = 0;
    for i = 1:viewNum
        tmp = (X{i} - U{i}*V')*W{i};
        ff = ff + sum(sum(tmp.^2));
    end
    if abs((ff-f)/f)<1e-4 | abs(ff-f)>1e100 | time == 30	
        break;
    end
    f = ff; 
end
end


