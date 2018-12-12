function R=fRot(Axis,alpha);
    switch Axis
        case 'x'
            R=fRotx(alpha);
        case 'y'
            R=fRoty(alpha);
        case 'z'
            R=fRotz(alpha);
    end
