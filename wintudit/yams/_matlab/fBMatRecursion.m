function Bi=fBMatRecursion(Bp,Bhat,R0p,r_pi)
    % Recursive formulae for B' and Bhat
    np=size(Bp,2);
    ni=size(Bhat,2);

    Bi=zeros(6,ni+np);

    for j=1:np
        Bi(1:3,j) = Bp(1:3,j)+cross(Bp(4:6,j),r_pi);
        Bi(4:6,j) = Bp(4:6,j);
    end

    if ni>0
        Bi(1:3,np+1:end)=R0p*Bhat(1:3,:);
        Bi(4:6,np+1:end)=R0p*Bhat(4:6,:);
    end

end
