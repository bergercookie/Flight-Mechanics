function [rho,aspeed,temp,press] = stdatm(altg)
%C...................................................950915..............
%c     Given altitude, compute temperature, density, pressure,
%c     and speed of sound using the standard atmosphere
%C........|.........|.........|.........|.........|.........|.........|..
global g

radius=6356766.0;
r=287.0;

% Compute geopotential altitude (altg given in km)

alt=1000*altg*radius/(radius+altg*1000);

%     Gradient region

if alt <= 11000.0
        aslope=-6.5D-3;
        alt1=0;
        temp1=288.16D0;
        press1=1.01325D5;
        rho1=1.2250D0;
        g=9.80665D0;
        temp=temp1+aslope*(alt-alt1);
        press=press1*(temp/temp1)^(-g/(aslope*r));
        rho=rho1*temp1*press/(press1*temp);
end

%C     Isothermal region

if alt > 11000.0D0
        temp=216.66D0;
        alt1=11000.0D0;
        press1=0.22627D5;
        rho1=0.3638350D0;
        g=9.80665D0;
        press=press1*exp(-(g*(alt-alt1)/(r*temp)));
        rho=rho1*press/press1;
end

%C     We are not designing a space plane

      if alt > 25000.0D0
         disp('stop Altitude > 25 km')
      end

%C     The speed of sound

      aspeed=sqrt(1.4D0*r*temp);

