function [ML]=MicroLenses(SH,ML,WF)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This function calculates a microlenses array. It is supposed a square detector and no separation between microlenses.

INPUTS:
	SH = structure with the user defined configuration of the Shack-Hartman.
    WF = Incoming wavefront to the Shack-Hartmann
	
OUTPUTS
	ML.mask= binary mask which contains the microlenses mask.
	ML.radius= Radius of each microlense in pixels.
	ML.coor= Coordinates of the upper left corner of each valid microlense. A valid microlense is the one which is completely inside the whole pupil of the system.
	ML.erased = Coordinates of the upper left corner of each invalid microlense. For example, discarded microlenses because they are incomplete.
%}


ML=CreatePupilML(ML);% Here, the pupil of each microlense is created.
% the upper left corner of the square where each microlense is inscribed is defined .
secuencex=ML.surplus(1)+1:ML.spacing:SH.resolution-ML.spacing+1; 
[p,q] = meshgrid(secuencex, secuencex);
pairs = [p(:) q(:)];
Circles = zeros(SH.resolution, SH.resolution); % Aux Matrix
x11=zeros(1,length(pairs));
x22=x11;
y11=x11;
y22=x11;


% Creating the microlenses of lenstlet one by one:
for k = 1 : length(pairs)
        % find upper left corner:
        x1 = int16(pairs(k,1));
        x11(k)=x1;
        y1 = int16(pairs(k,2));
        y11(k)=y1;
        x2 = int16(x1 + ML.spacing - 1);
        x22(k)=x2;
        y2 = int16(y1 + ML.spacing - 1);
        y22(k)=y2;
        % Adding the tiny pupil and the phase of each
        Circles (y1:y2, x1:x2)= Circles(y1:y2,x1:x2) + ML.Pupil; % Aux Matrix
end


% Multiplication of the main mask and the microlenses mask. 
ML.AmplitudeMask = Circles.*SH.pupil;

% Now, only those microlenses which are completely inside the main pupil are selected: 
ML.coor=zeros(1,4);
ML.erased=zeros(length(y11),1);
for i=1:length(x11)
    if sum(sum(ML.AmplitudeMask(y11(i):y22(i), x11(i):x22(i))))<sum(sum(ML.Pupil))
        ML.AmplitudeMask(y11(i):y22(i), x11(i):x22(i))=0;
    else 
        ML.coor(end+1,1:4)=[y11(i),y22(i),x11(i),x22(i)];
        ML.erased(i)=1;
    end
end

ML.coor(1,:)=[]; %In the way these vectors are constructed, the first line is always zero.
ML.Radius_mm=SH.PixelSize*ML.radiusPixels;%Radius in mmm
ML.AmplitudeMask4Paint=ML.AmplitudeMask;
ML.AmplitudeMask4Paint(ML.AmplitudeMask4Paint==0)=NaN;

%Painting...
if SH.paint==1
    subplot(1,3,2)
    imshow(ML.AmplitudeMask,[]);
    title('Microlenses mask');
    set(gcf,'color','w');
    xlabel('pixels')
    ylabel('pixels')
    colorbar

    subplot(1,3,3)
    imshow(ML.AmplitudeMask4Paint.*WF,[])
    title('Incoming wavefront through lenslet (phase of microlenses is not yet applied)')
    set(gcf,'color','w');
    xlabel('pixels')
    ylabel('pixels')
    drawnow();
    colorbar
end
