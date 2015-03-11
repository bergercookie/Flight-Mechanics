function [xdot]=fplmod(tval,xstate)
% This version with all variables as states
% Use a wrapper for most types of analysis
% 2007-05-15: alfdot and tetdot statements corrected
% 2008-01-07: Draken based on VERY preliminar longitudinal data
% 2008-01-10: New mass data from Bo Nilsson
% 2008-01-17: Improved low speed data for -50 < alfa < +140 deg
% 2008-01-18: Reduced to longitudinal DOFs only
% 2009-02-09: Some simplifications


% globals
global xcg

% Some constants in use

zero=0.0;half=0.5;one=1.0;two=2.0;
r2d=57.29577951308232;kappa=1.4;
expg1=(kappa-one)/kappa;
expg2=kappa/(kappa-one);
aspeed0=340.269;press0=101325.0;
tusen=1000.0;

% Aircraft data

      xref=10.0;             % Reference point for aerodynamics data (m)
      yref=0.0;
      zref=0;
      cbar = 5.3191;         % Reference chord (m)
      span =    9.4;         % Span (m)
      Sref =   50.0;         % Reference area (m*m)
      gacc =   9.81;         % Acceleration of gravity (m/(s*s))
      tepsr=-5.0*pi/180;     % Thrust installation angle in radians (5 deg)

%     Aircraft mass without fuel

      amass0 = 8180;

%     Define all the variables

      u        = xstate(1);    % m/s
      w        = xstate(2);    % m/s
      qrate    = xstate(3);    % rad/s
      theta    = xstate(4);    % rad
      xdiste   = xstate(5);    % m
      altm     = xstate(6);    % m
      fuelmass = xstate(7);    % kg
      deltae   = xstate(8);    % rad
      deltap   = xstate(9);    % 0 (flight idle) to 1 (full thrust) 

%     Aircraft mass

      amass = amass0+fuelmass;

%     Transformation matrix body to earth

      sintet=sin(theta);
      costet=cos(theta);

%     Compute wind axis state

      airspd = sqrt( u*u + w*w );
      alfa   = atan(w/u);

%     More extra variables

      sinalf=sin(alfa);
      cosalf=cos(alfa);

%     These should be a function of the fuel level

%       xcg=10.1;
      zcg=0;
      z_eng=0;                 % Engine c.g. offset

% Gunnar Strangs rapport

      aIy    = 75300.0;       % Pitching moment of inertia

%     Define some temporary variables

      altkm=altm/1000;       % Altitude in km

%     Get atmospheric data

      [rho,aspeed,temp,press] = stdatm(altm);

%     Compute Mach number

      rmach=airspd/aspeed;
      rmach2=rmach*rmach;
      qdyn=rho*airspd*airspd/two;

%     Compute calibrated (indicated) air speed
%     First, compute the total pressure that the aircraft should sense

      ptot=press*(one+((kappa-one)/two)*rmach2)^expg2;
      qc=ptot-press;
      rmtem=sqrt( (two/(kappa-one))*( (qc/press0 + one)^expg1 - one ) );
 
%     Indicated (calibrated) airspeed and derivatives

      vcal=rmtem*aspeed0;

%     Nondimensional pitch rate

      qbar=qrate*cbar/(2*airspd);

%     Engine data (Simple scaling of full thrust no afterburner)

      iebk=0;

      [thrust,fuelb]=rm6cal(rmach,altkm,iebk);

      thrust=thrust*deltap;
      fuelb=fuelb*deltap;

%     New aerodata, low speed only but up to high alfa

      [CD,CL,Cm]=clongm02(alfa,deltae);
      CLq=2.7;
      CL=CL+CLq*qbar;
      rlift = qdyn*Sref*CL;
      drag = qdyn * Sref *CD;
      Cmq=-1.25;
      Cm = Cm + Cmq*qbar;

      Cmadot = -0.25;

%     Stevens and Lewis model as on p. 123-124.
%     Forces in body axis system

      xforce = thrust - drag*cosalf  + rlift*sinalf;
      zforce =        - rlift*cosalf - drag*sinalf;

%     Force equation(s) body axis

      udot = - qrate*w + xforce/amass - gacc*sintet;
      wdot = qrate*u + zforce/amass + gacc*costet;

%     Time derivative of airspeed

      aspdot = (u*udot + w*wdot)/airspd;

%     Time derivative of angle-of-attack

      alfdot = (cosalf*cosalf)*(wdot*u-w*udot)/(u*u);

%     Moment equation(s)

pitmom = (qdyn*Sref*cbar*(Cm + Cmadot*alfdot*cbar/(2*airspd))+thrust*z_eng);

%******************************%
%     Adjust for c.g. position
%******************************%
      pitmom = pitmom + rlift* ( cosalf*(xcg-xref) - sinalf*(zcg-zref) ) ...
                      + drag*  ( cosalf*(zcg-zref) + sinalf*(xcg-xref) );

%     Vehicle symmetry is assumed

      qdot = ( pitmom )/aIy;

%     Where we are equations

      xedot  =  costet*u + sintet*w;
      altdot = - ( -sintet*u + costet*w);

      tetdot = qrate ;

%     State equations

      xdot(1:13)=0;

      xdot(1) = udot;
      xdot(2) = wdot;
      xdot(3) = qdot;
      xdot(4) = tetdot;
      xdot(5) = xedot;
      xdot(6) = altdot;
      xdot(7) = -fuelb;

%     Some useful additional information

%      nz=qdyn*Sref*CL/(amass*gacc);
      nz=-zforce/(amass*gacc);

%     Extra output variables, add more if you like
% I 'll have to modify this !
% add more variables
      xdot(8)=airspd;
      xdot(9)=alfa;
      xdot(10)=vcal;
      xdot(11)=rmach;
      xdot(12)=CL;
      xdot(13)=nz;
