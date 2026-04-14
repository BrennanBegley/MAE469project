function [a, eVec, eMag, inc, RAAN, argPeri, trueAnom] = orbitalelementscalc(r_ijk, v_ijk, mu)
%ORBITALELEMENTSCALC Compute classical orbital elements from inertial state.
%   [a, eVec, eMag, inc, RAAN, argPeri, trueAnom] = orbitalelementscalc(r,v,mu)
%   returns the classical orbital elements associated with the state vector.
%
% Outputs:
%   a         Semi-major axis
%   eVec      Eccentricity vector
%   eMag      Eccentricity magnitude
%   inc       Inclination [rad]
%   RAAN      Right ascension of ascending node [rad]
%   argPeri   Argument of periapsis [rad]
%   trueAnom  True anomaly [rad]

    rMag = norm(r_ijk);
    vMag = norm(v_ijk);

    hVec = cross(r_ijk, v_ijk);
    hMag = norm(hVec);

    kHat = [0 0 1];
    nVec = cross(kHat, hVec);
    nMag = norm(nVec);

    eVec = (1/mu) * ((vMag^2 - mu/rMag) * r_ijk - dot(r_ijk, v_ijk) * v_ijk);
    eMag = norm(eVec);

    specificEnergy = vMag^2/2 - mu/rMag;
    a = -mu / (2 * specificEnergy);

    inc = acos(dot(kHat, hVec) / hMag);

    if nMag == 0
        RAAN = NaN;
    else
        RAAN = acos(dot([1 0 0], nVec) / nMag);
        if nVec(2) < 0
            RAAN = 2*pi - RAAN;
        end
    end

    if eMag == 0 || nMag == 0
        argPeri = NaN;
    else
        argPeri = acos(dot(nVec, eVec) / (nMag * eMag));
        if eVec(3) < 0
            argPeri = 2*pi - argPeri;
        end
    end

    if eMag == 0
        trueAnom = acos(dot([1 0 0], r_ijk) / rMag);
    else
        trueAnom = acos(dot(eVec, r_ijk) / (eMag * rMag));
        if dot(r_ijk, v_ijk) < 0
            trueAnom = 2*pi - trueAnom;
        end
    end
end
