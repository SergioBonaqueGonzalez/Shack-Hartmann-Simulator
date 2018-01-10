function ML=CreatePupilML(ML)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergiob@wooptix.com
This program calculates a variery of geometry for pupils of the microlenses
array of a Shack-Hartmann Sensor
%}
xp=linspace(-1,1,ML.spacing);
SmallCircle = ones(ML.spacing, ML.spacing);
[X,Y]=meshgrid (xp,xp);
[rho]=sqrt(X.^2+Y.^2);
[a,b]=size(rho);
for i=(1:a);
    for j=(1:b);
        if rho(i,j) > 1 ;
            SmallCircle(i,j)=0;
        end;
    end;
end;

if ML.Prop==0
    ML.Pupil=SmallCircle;
    ML.Phase=ones(length(SmallCircle));
else
    if ML.Geometry==1 
        ML.Pupil=SmallCircle; %Pupil
    else
        ML.Pupil=ones(length(SmallCircle));
    end
end

