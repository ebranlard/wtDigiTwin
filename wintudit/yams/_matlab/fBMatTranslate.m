function Bi=fBMatTranslate(Bp,r_pi)
    % Translate a B matrix
    Bi=zeros(size(Bp));
    for j=1:size(Bp,2)
        Bi(1:3,j) = Bp(1:3,j)+cross(Bp(4:6,j),r_pi);
        Bi(4:6,j) = Bp(4:6,j);
    end
end
