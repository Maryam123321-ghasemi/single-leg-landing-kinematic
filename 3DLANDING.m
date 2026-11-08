clc;
clear all;
%% kinematics
%load files
filename=('L1-2.txt');
dat=importdata(filename);           %load trial matrix
data_L=dat.data(3:end,:);             % drop landing

filename1=('ana-1.txt');
dat1=importdata(filename1);
data_ana=dat1.data(3,:);
%% Filtering
fs=200;
fc=7;
order=4;
[b,a]=butter(order,fc/(fs/2),'low'); 
data1=filtfilt(b,a,data_L);
%time
data1(:,1)=data_L(:,1);

%define Anatomic points
abd=data_ana(:,2:4);
r_a=data_ana(:,5:7);            %anatomic asis L&R
l_a=data_ana(:,8:10); 
l5=data_ana(:,11:13);           %l5 ana (s1)
ca_1=data_ana(:,14:16);         %cluster anat c1-c6  
ca_2=data_ana(:,17:19);
ca_3=data_ana(:,20:22);
kl=data_ana(:,23:25);            %Lateral knee 
km=data_ana(:,26:28);            %Medial knee 
ca_4=data_ana(:,29:31);
ca_5=data_ana(:,32:34);
ca_6=data_ana(:,35:37);
ml=data_ana(:,38:40);      %lateral malleolus
mm=data_ana(:,41:43);       %medial malleolus
h=data_ana(:,44:46);       %heel
m2=data_ana(:,47:49);
m5=data_ana(:,50:52);      %met 5
to=data_ana(:,53:55);       %toe

%define Tracking points
Asis_r=data1(:,2:4);
Asis_l=data1(:,5:7);
L5=data1(:,8:10);
c1=data1(:,11:13);          %markers clusters
c2=data1(:,14:16);
c3=data1(:,17:19);
c4=data1(:,20:22);
c5=data1(:,23:25);
c6=data1(:,26:28);
he=data1(:,29:31);        %heel marker
met_2=data1(:,32:34);     %met2 marker
met_5=data1(:,35:37);     %met5 marker

%% Joints centers & Segments CM
o_p=(r_a+l_a)/2;   %pelvis anatomic coordinate system
oo=l5;
x_p=r_a-o_p;                 
yy_p=oo-o_p;
z_p=cross(yy_p,x_p);
y_p=cross(z_p,x_p);

R=norm(r_a-l_a);                   % |asis R- asis L |            
Pr_h=[0.36*R -0.19*R -0.30*R];  %right hip formula in pelvis anatomic coordinate system.  
Pl_h=[-0.36*R -0.19*R -0.30*R];  %left hip formula in pelvis anatomic coordinate system.
rp_h=[x_p(:)/norm(x_p(:)), y_p(:)/norm(y_p(:)),z_p(:)/norm(z_p(:))];%rotational matrix pelvis

hip_R=(Pr_h*rp_h)+o_p;            %hip_r joint center   
hip_L=(Pl_h*rp_h)+o_p;            %hip_l joint center  
hips=(hip_R+hip_L)/2;             %hip joint center    

HIP=hip_R;                                 
k=(kl+km)/2;                   % anatomic knee joint center
A=(ml+mm)/2;                  % anatomic ankle joint center=maleolus l+maleolus m/2

%men
cmpelvis=0.3885*abd+0.6115*hips;        %pelvis cm 
cmthigh=0.5905*HIP+0.4095*k;             %thigh cm                
cmshank=0.5605*k+0.4395*A;               %shank cm                
cmfoot=0.5585*h+0.4415*to;                %foot cm                                             

%% define anatomical axes
 %thigh
z_th=HIP-k;
x_thth=km-k;
y_th=cross(x_thth,z_th);
x_th=cross(y_th,z_th);
%shank        
z_sh=k-A;                      
x_shsh=mm-A;
y_sh=cross(x_shsh,z_sh);                       
x_sh=cross(y_sh,z_sh);
%foot 
y_fo=to-h;                         %y=toe-heel               
x_fofo=m5-h;
z_fo=cross(x_fofo,y_fo);                       
x_fo=cross(y_fo,z_fo); 
% norm anatomical axes
r_P=[x_p/norm(x_p);y_p/norm(y_p);z_p/norm(z_p)];         
r_T=[x_th/norm(x_th);y_th/norm(y_th);z_th/norm(z_th)];  
r_S=[x_sh/norm(x_sh);y_sh/norm(y_sh);z_sh/norm(z_sh)];
r_F=[x_fo/norm(x_fo);y_fo/norm(y_fo);z_fo/norm(z_fo)];

%% Transform tracking Markers from GCS to ACS 
%pelvis
O=o_p-o_p;                   
anat_O=O*r_P;              

R_A=r_a-o_p;             %defined at line 54.x_p.
%R_A=r_a-cmpelvis;      
anat_r_a=R_A*r_P;

L_A=l_a-o_p;
%L_A=l_a-cmpelvis;
anat_l_a=L_A*r_P;

L_5=l5-o_p;
%L_5=l5-cmpelvis;
anat_l5=L_5*r_P;          

 %thigh anatomic clusters
C1=ca_1-k;          
%C1=ca_1-cmthigh;
anat_c1=C1*r_T;

C2=ca_2-k;
%C2=ca_2-cmthigh;
anat_C2=C2*r_T;

C3=ca_3-k;
%C3=ca_3-cmthigh;
anat_c3=C3*r_T;

%shank anatomic clusters
C4=ca_4-A;              
%C4=ca_4-cmshank;
anat_c4=C4*r_S;

C5=ca_5-A;          
%C5=ca_5-cmshank;
anat_C5=C5*r_S;

C6=ca_6-A;
%C6=ca_6-cmshank;     
anat_c6=C6*r_S;

%foot
H=h-h;
%H=h-cmfoot;
anat_h=H*r_F;

M2=m2-h;
%M2=m2-cmfoot;
anat_met2=M2*r_F;

M5=m5-h;
%M5=m5-cmfoot;
anat_met5=M5*r_F;      


%% Define Constant Matrices(definition of MCS's at ACS)

%pelvis
YY_A_P=anat_l5-anat_l_a;       %Y'
X_A_P=anat_r_a-anat_l_a;
Z_A_P=cross(YY_A_P,X_A_P);
Y_A_P=cross(Z_A_P,X_A_P);

%thigh
Z_A_T=anat_c1-anat_C2;
YY_A_T=anat_c3-anat_C2;
X_A_T=cross(YY_A_T,Z_A_T);
Y_A_T=cross(Z_A_T,X_A_T);

%shank
Z_A_S=anat_c4-anat_C5;
YY_A_S=anat_c6-anat_C5;
X_A_S=cross(YY_A_S,Z_A_S);
Y_A_S=cross(Z_A_S,X_A_S);

%foot
XX_A_F=anat_met5-anat_h;            %X'_ANAT_M5  
Y_A_F=anat_met2-anat_h;
Z_A_F=cross(XX_A_F,Y_A_F);
X_A_F=cross(Y_A_F,Z_A_F);

% norm marker axes
C_R_P=[X_A_P/norm(X_A_P);Y_A_P/norm(Y_A_P);Z_A_P/norm(Z_A_P)];%pelvis

C_R_T=[X_A_T/norm(X_A_T);Y_A_T/norm(Y_A_T);Z_A_T/norm(Z_A_T)];%thigh

C_R_S=[X_A_S/norm(X_A_S);Y_A_S/norm(Y_A_S);Z_A_S/norm(Z_A_S)];%shank

C_R_F=[X_A_F/norm(X_A_F);Y_A_F/norm(Y_A_F);Z_A_F/norm(Z_A_F)];%FOOT

%% Define marker trial axes(from tracking)

%pelvis
x_pe=Asis_r-Asis_l;         
y_pp=L5-Asis_l;                      
z_pe=cross(y_pp,x_pe);
y_pe=cross(z_pe,x_pe);

%thigh vec
z_t=c1-c2;                 %cluster 1:3%marker coordinate system
x_tt=c3-c2;
y_t=cross(z_t,x_tt);
x_t=cross(y_t,z_t);

%shank vec
z_s=c4-c5;                 %cluster 4:6%marker coordinate system
x_ss=c6-c5;
y_s=cross(z_s,x_ss);                       
x_s=cross(y_s,z_s);

%foot vec
y_f=met_2-he;           %heel, met 5&2. marker coordinate system           
x_ff=met_5-he;
z_f=cross(x_ff,y_f);
x_f=cross(y_f,z_f);


%pelvis
for i=1:length(z_pe)
r_pe(3*i-2:3*i,:)=[x_pe(i,:);y_pe(i,:);z_pe(i,:)]; 
end

 %thigh 
for i=1:length(z_t)
r_th(3*i-2:3*i,:)=[x_t(i,:);y_t(i,:);z_t(i,:)];
%%r_th1(:,:,i)=[x_t(i,:)/norm(x_t(i,:)); y_t(i,:)/norm(y_t(i,:)); z_t(i,:)/norm(z_t(i,:))];  
end

 %shank
for i=1:length(z_s)
    r_sh(3*i-2:3*i,:)=[x_s(i,:);y_s(i,:);z_s(i,:)];
   %% r_sh1(:,:,i)=[x_s(i,:)/norm(x_s(i,:)); y_s(i,:)/norm(y_s(i,:)); z_s(i,:)/norm(z_s(i,:))];  
end

%foot 
for i=1:length(z_f)
    r_foo(3*i-2:3*i,:)=[x_f(i,:);y_f(i,:);z_f(i,:)];  
end
%% Rotational Matrices
%every segment's rotational matrix in global
for i=1:length(z_pe)          %pelvis 
  
  R_PELVIS(3*i-2:3*i,:)=C_R_P*r_pe(3*i-2:3*i,:);               
end

for i=1:length(z_t)          %thigh
 R_THIGH(3*i-2:3*i,:)=C_R_T*r_th(3*i-2:3*i,:);
end

for i=1:length(z_s)          %shank
  R_SHANK(3*i-2:3*i,:)=C_R_S*r_sh(3*i-2:3*i,:);
end

for i=1:length(z_f)          %foot
 R_FOOT(3*i-2:3*i,:)=C_R_F*r_foo(3*i-2:3*i,:);  
end  
%% angles calculations

%pelvis
for i=1:length(z_pe);
alpha_pe(i,:)=atan2d(-R_PELVIS(3*i,2),R_PELVIS(3*i,3));
beta_pe(i,:)=atan2d(R_PELVIS(3*i,1),sqrt(R_PELVIS(3*i-2,1)^2+R_PELVIS(3*i-1,1)^2));
gamma_pe(i,:)=atan2d(-R_PELVIS(3*i-1,1),R_PELVIS(3*i-2,1));
end
%thigh
for i=1:length(z_t);
alpha_t(i,:)=atan2d(-R_THIGH(3*i,2),R_THIGH(3*i,3));
beta_t(i,:)=atan2d(R_THIGH(3*i,1),sqrt(R_THIGH(3*i-2,1)^2+R_THIGH(3*i-1,1)^2));
gamma_t(i,:)=atan2d(-R_THIGH(3*i-1,1),R_THIGH(3*i-2,1));
end
  %shank
for i=1:length(z_s);
alpha_s(i,:)=atan2d(-R_SHANK(3*i,2),R_SHANK(3*i,3));
beta_s(i,:)=atan2d(R_SHANK(3*i,1),sqrt(R_SHANK(3*i-2,1)^2+R_SHANK(3*i-1,1)^2));
gamma_s(i,:)=atan2d(-R_SHANK(3*i-1,1),R_SHANK(3*i-2,1));
end

for i=1:length(z_f);
alpha_f(i,:)=atan2d(-R_FOOT(3*i,2),R_FOOT(3*i,3));
beta_f(i,:)=atan2d(R_FOOT(3*i,1),sqrt(R_FOOT(3*i-2,1)^2+R_FOOT(3*i-1,1)^2));
gamma_f(i,:)=atan2d(-R_FOOT(3*i-1,1),R_FOOT(3*i-2,1));
end