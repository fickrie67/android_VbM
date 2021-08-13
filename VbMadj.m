clc, clear all
format long
file1 = 'marker9.xlsx';
file2 = 'vslam_pose9.xlsx';
%file3 = 'marker9.xlsx';

A = xlsread(file1);
B = xlsread(file2);
C = xlsread(file1, 'marker9');
D = xlsread(file2, 'vslam_pose9');
E = xlsread(file1, 'Sheet1');
F = xlsread(file2, 'Sheet1');

%%%%%%%%%%%%%%%%%%%
%%adj%%
%marker adj
 xm=A(:,15);
 ym=A(:,16);
 zm1=A(:,17);
 zm=(zm1+0.067)-0.183;
 t=A(:,26);
 
%ptam adj
 xv=B(:,7) ;
 yv=B(:,8);
 zv=B(:,9);
 
 
 time1=B(:,1);
 t1=uint64(1625714747997840000);
 %time=uint64(t1);
 %time2=B(:,2);
 date = datetime(t1,'ConvertFrom','epochtime','TicksPerSecond',1e9,'Format','dd-MMM-yyyy HH:mm:ss.SSSSSSSSS');
%date_time = datestr(unix_time_pose./86400 + datenum(1970,1,1));
 t2 = uint64(1545390864126080000);
d = datetime(t2,'ConvertFrom','epochtime','TicksPerSecond',1e9,'Format','dd-MMM-yyyy HH:mm:ss.SSSSSSSSS');
 
 o=A(:,24);
 l=A(:,25);
%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%
%marker
 x2=C(:,15);
 y2=C(:,16);
 z2a=C(:,17);
 z2=(z2a+0.067)-0.183;
 
 %ptam
 x1=D(:,7) ;
 y1=D(:,8);
 zl=D(:,9);
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%plot error
%marker
 xmarker=E(:,15);
 ymarker=E(:,16);
 zmarker1=E(:,17);
 zmarker=(zmarker1+0.067)-0.183;
 timeA=E(:,1);
 timeB = timeA(~isnan(timeA));
 timeC=uint64(timeB);
 dateA = datetime(timeC,'ConvertFrom','epochtime','TicksPerSecond',1e9,'Format','HH:mm:ss');
 
 timeD = E(:,27);
 timeE = timeD(~isnan(timeD));
 
 %ptam
 xvslam=F(:,7) ;
 yvslam=F(:,8);
 zvslam=F(:,9);
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%design matrix%%%%%%%
i=1:160;
AA=[xv(i) o(i) o(i) l(i) o(i) o(i)
    o(i) yv(i) o(i) o(i) l(i) o(i)
    o(i) o(i) zv(i) o(i) o(i) l(i)];
BB=[xm(i) ;ym(i); zm(i)];

%disp (AA)
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%calculate parameters%%
%X=(AA'*AA)^-1*(AA'*BB);
 NN=(AA'*AA);
 N=NN^-1;
 n=AA'*BB;
 X=N*n
 
 %%%%%parameters%%%%%%%%%
 sx=X(1,1);
 sy=X(2,1);
 sz=X(3,1);
 tx=X(4,1);
 ty=X(5,1);
 tz=X(6,1);
 
 
 %%%calculate vslam pos%%%
 xma=xm*100;
 yma=ym*100;
 zma=zm*100;
 
 xxa=(xv*sx+tx)*100;
 yya=(yv*sy+ty)*100;
 zza=(zv*sz+tz)*100;
 
 %%%%%%%%%%%%%%%%%%%%%%%%
 xmb=x2*100;
 ymb=y2*100;
 zmb=z2*100;
 
 xx=(x1*sx+tx)*100;
 yy=(y1*sy+ty)*100;
 zz=(zl*sz+tz)*100;
 
 %%%%%%%%%%%%%%%%%%%%%%%%
 xmc=xmarker*100;
 ymc=ymarker*100;
 zmc=zmarker*100;
 
 xxc=(xvslam*sx+tx)*100;
 yyc=(yvslam*sy+ty)*100;
 zzc=(zvslam*sz+tz)*100;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%rmse%%%%%%%%%%%%%%
xvar1=xmc-xxc;
xvar=xvar1(~isnan(xvar1));
yvar1=ymc-yyc;
yvar=yvar1(~isnan(yvar1));
zvar1=zmc-zzc;
zvar=zvar1(~isnan(zvar1));

 vxmean=(mean(xvar));
 vymean=(mean(yvar));
 vzmean=(mean(zvar));
 
 xqr=xvar.^2;
 yqr=yvar.^2;
 zqr=zvar.^2;
 radius=sqrt(xqr+yqr+zqr);
 
 Xmax = max(xvar(:));
 Ymax = max(yvar(:));
 Zmax = max(zvar(:));
 Rmax = max(radius(:));
 

 Xmin = min(xvar(:));
 Ymin = min(yvar(:));
 Zmin = min(zvar(:));
 Rmin = min(radius(:));
  
 meanX=mean(xvar);
 meanY=mean(yvar);
 meanZ=mean(zvar);
 meanR = mean(radius(:));
 
 Rselisih=radius-meanR;
 Rabs=abs(Rselisih);
 Rsum= sum(Rabs);
 Rabsmean=Rsum/446;
 
 stdX = std(xvar(:));
 stdY = std(yvar(:));
 stdZ = std(zvar(:));
 stdR = std(radius(:));
 
 meanxqr=mean(xqr);
 meanyqr=mean(yqr);
 meanzqr=mean(zqr);

 Xrmse=sqrt(meanxqr);
 Yrmse=sqrt(meanyqr);
 Zrmse=sqrt(meanzqr);
 %Rrmse=sqrt(
 
 
 
 RMSE = [Xrmse Yrmse Zrmse];
 
 Vv=[(xvar)/100 ;(yvar)/100 ;(zvar)/100];
 
 sigv=(Vv'*Vv)/(160)
 Cxv=sigv*N;
 SVall=sqrt(diag(Cxv))
%%%%%parameters%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%plot%%%%%%%%%%%%%%
 
 figure(1)
 plot(timeE,xvar, 'LineWidth',2);
 grid on
 %title ('Position Drift of VbM against Groundtruth')
 hold on
 plot(timeE,yvar, 'LineWidth',2);
 plot(timeE,zvar, 'LineWidth',2);
 plot(timeE,radius, 'LineWidth',2);

 xlabel('Time [s]')
 %xticks([0 25 50 75 100 125 150])
 ylabel('Drift Error [cm]')
 %yticks([-12 -8 -4 0 4 8 12])
 
 ax = gca;
 ax.FontSize = 13;
legend('X','Y','Z','R');

figure(2)
plot3(xmb,ymb,zmb,'-',xx,yy,zz,'-','LineWidth',3);
grid on
%title ('VbM Trajectories Underwater')
%xticks([-10 -6 -2 2 6 10])
%yticks([-10 -8 -6 -4 -2 0 ])
%zticks([-65 -60 -55])
zlim([-5 0])
grid on

xlabel('X [cm]')
ylabel('Y [cm]')
zlabel('Z [cm]')
ax = gca;
 ax.FontSize = 14;

legend('Ground-truth','VbM');
%%%%%%%%%%%%%%%%%%%%%%%%%