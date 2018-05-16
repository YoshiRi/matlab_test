function [hl,hr] = ImageObs(f,X,Z,Cta,W,H)

hl = f*H/(Z-W/2*cos(Cta));
hr = f*H/(Z+W/2*cos(Cta));

end