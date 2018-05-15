clc
clear
%In this code we try to solve unsteady heat transfer problem over
%rectangular domain Toward the steady Soltion.
%Right and Left Walls are kept at T0, T1 Top and bottom walls are insulated
%Left and  Right Boundries are kept in T0 and T1, respectively
%Top and Buttom Boundries are Zero Flux
%BCS Can be change in the Bcs Function
%By Mohammad Aghakhani, 2012

%  *************************Top******************************
%  *........................................................*
%  *........................................................*
% Left--->T0......................................T1<---Right
%  *........................................................*
%  *........................................................*
%  ***********************Bottom*****************************

%-----------------Inputs-----------
L=2;   %Channel Lenght
H=1;   %Chnnel with
T0=300; %Left wall Temprature
T1=50;%Right wall temperature
t0=100;   %iniial Values
alpha=0.00023; %Thermal diffusivity
m=100;  % No. of points along Top & Bottom
n=100 ; %No. of point along Left & Right sides
IT=1;    %Current iteration No.
MIT=100000; %Maximum allowabe iteration
Dt=0.22; %time step
eps=0.51e-3; %error
errT=1000; %Error in two con. time ste p
intAnimT=100; %interval of Animation
%------------------------------------

%Coordinate of nodes
X=zeros(n,m);
Y=zeros(n,m);

%Solution varibles @ Current  iteration
T=zeros(n,m);

%Solution varibles @ Previous iteration
Told=zeros(n,m);
%Animamtion congfig
%Set Figures for Animation
scrsz = get(0,'ScreenSize');
%Contour
CT=figure('Name','Animation of the Temperature Contour','NumberTitle','off','OuterPosition',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
vidObjC = VideoWriter('TemperatureContour.avi');
vidObjC.Quality = 100;
vidObjC.FrameRate = 5;
open(vidObjC);
%Profie
PT=figure('Name','Animation of the Temperature Profile','NumberTitle','off','OuterPosition',[scrsz(1) scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
vidObjP = VideoWriter('TemperatureProfile.avi');
vidObjP.Quality = 100;
vidObjP.FrameRate = 5;
open(vidObjP);
%Error-Convergence
RT=figure('Name','Animation of the Convrgence','NumberTitle','off','OuterPosition',[scrsz(1) scrsz(2) scrsz(3)/2 scrsz(4)/2]);
vidObjR = VideoWriter('ConvergenceHistory.avi');
vidObjR.Quality = 100;
vidObjR.FrameRate = 5;
open(vidObjR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grid Genration
[X,Y,dL,dH]=Grid(m,n,L,H);
F0=alpha*Dt/(dL*dH);
fprintf(1,'Fourier Number  =  %2.5f\n',F0);
%%disp('Press any Key to continue');
% pause
if(F0>=.25)
    disp('Fourier Number is too High, Please decrease Time Step');
    return;
end

%Initiate the solution
[T]=initiate(n,m,t0);

%Set boundry condition for first time
[T]=Bcs(n,m,T,T0,T1);
Told=T;
error=zeros(1,MIT);
%Begin Iteration
while((IT<MIT)&&(errT>eps))
    %shift solution from old iteration
    Told=T;
    [T]=Bcs(n,m,T,T0,T1);
    %Caluate Right hand side(aplha*D^T/DX^2+*D^T/DY^2)
    %RHST=alpha*(DDXX(n,m,dL,T)+DDYY(n,m,dH,T));
    %Explicit Euler Method
    for i=2:n-1
        for j=2:m-1
            DTDXX=(  T(i,j+1)-2.*T(i,j)+T(i,j-1)  )/(dL.*dL);
            DTDYY=(  T(i+1,j) -2*T(i,j)+T(i-1,j)  )/(dH.*dH);
            T(i,j)=T(i,j)+alpha *( DTDXX+DTDYY)*Dt;
        end
    end
    %Calculate Errors
    errT=max(max(abs((T-Told))));
    error(IT)=errT;
    fprintf(1,'IT=%i  Time=%2.6e  Error=%2.6e\n',IT,IT*Dt,errT);
    %animatation
    if mod(IT,intAnimT)==0 || IT==1
        % Animation
        %----Plot  T Contour
        figure(CT);
        [C,h]=contourf(X,Y,T);
        xlabel('X Coordinate');
        ylabel('Y Coordinate');
        title(strcat('Temperature Contour  ','  (','Time Step=',num2str(IT),' ,Time= ',num2str(Dt*(IT)),')'));
        clabel(C,h)
        drawnow
        writeVideo(vidObjC, getframe(gca));
        %Profile Animation
        figure(PT);
        %----Plot  T Profile
        plot(X(floor(n/2),:),T(floor(n/2),:))
        xlabel('X Coordinate');
        ylabel('Y Coordinate');
        title(strcat('Temperature Profile at Mid-Secton ','  (','Time Step=',num2str(IT),' ,Time= ',num2str(Dt*(IT)),')'));
        clabel(C,h)
        drawnow
        writeVideo(vidObjP, getframe(gca));
        %Convergence Animation
        figure(RT);
        %----Plot  T Contour
        semilogy(1:IT,error(1:IT),'- r');
        xlabel('Iteration');
        ylabel('Error');
        title('Convergence History');
        drawnow
        writeVideo(vidObjR, getframe(gca));

        pause(0.5);
    end
    IT=IT+1;
end
if(errT<eps)
    fprintf(1,'Converged in %i Iterations',IT);
else
    disp('Maximum Iteration Number Reached');
    semilog(1:IT,error(1:IT),'- r');
    xlabel('Iteration');
    ylabel('Error');
    title('Convergence History');
    return;
end
%Postprocces
semilogy(1:IT,error(1:IT),'- r','LineWidth',2);
xlabel('Iteration');
ylabel('Error');
title('Convergence History');
% Countour of Temperature
figure
[C1,h1] = contourf(T,30);
text_handle = clabel(C1,h1,'manual');
colorbar
xlabel('x')
ylabel('y')
title(strcat('Temperature Contour @',num2str(IT*Dt)));
drawnow
% Surface of Temperature
figure
surf(X,Y,T);
title(strcat('Temperature Surface @',num2str(IT*Dt)));
xlabel('x')
ylabel('y')
zlabel('Temperature')
axis fill;
% Profile of Temperature at mid-Section
figure
plot(X(floor(n/2),:),T(floor(n/2),:),'LineWidth',2)
%Animation Output
close(vidObjC);
close(vidObjp);
close(vidObjR);
% winopen('TemperatureContour.avi')
display('AVI Movie(s) for Temperature Contour exported to the Current Directory')

