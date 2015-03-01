function [cd,cl,cm]=clongm02(alfa,deltae)
% Nondimensional drag coefficient for low Machnumbers (M=0.2)
% Diagram 20-23
% alfa and deltae in radians

global cdtab cltab cmtab

alfadeg=max(-50,min(140,alfa*180/pi));

dedeg=deltae*180/pi;

cd_de0=interp1(cdtab(:,1),cdtab(:,3),alfadeg);
cl_de0=interp1(cltab(:,1),cltab(:,3),alfadeg);
cm_de0=interp1(cmtab(:,1),cmtab(:,3),alfadeg);

if dedeg<0
  cd_de25=interp1(cdtab(:,1),cdtab(:,2),alfadeg);
  cd=cd_de0*(1+dedeg/25) + cd_de25*(-dedeg/25);

  cl_de25=interp1(cltab(:,1),cltab(:,2),alfadeg);
  cl=cl_de0*(1+dedeg/25) + cl_de25*(-dedeg/25);

  cm_de25=interp1(cmtab(:,1),cmtab(:,2),alfadeg);
  cm=cm_de0*(1+dedeg/25) + cm_de25*(-dedeg/25);
else
  cd_de10=interp1(cdtab(:,1),cdtab(:,4),alfadeg);
  cd=cd_de0*(1-dedeg/10) + cd_de10*(dedeg/10);

  cl_de10=interp1(cltab(:,1),cltab(:,4),alfadeg);
  cl=cl_de0*(1-dedeg/10) + cl_de10*(dedeg/10);

  cm_de10=interp1(cmtab(:,1),cmtab(:,4),alfadeg);
  cm=cm_de0*(1-dedeg/10) + cm_de10*(dedeg/10);
end

