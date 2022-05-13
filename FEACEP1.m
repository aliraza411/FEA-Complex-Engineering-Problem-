clc; clear all;
EE=210e9;                       %Modulus of elasticity
AA=80e-3;                       %Area of cross section for all elements m2
II=1.2e-4;                      %Area moment of inertia for all elements m4
LL1=3;                          %Length of element 1 m
LL2=6;                          %Length of element 2 m
Theta_1=0;                      %Angle in degree
Theta_2=90;                     %Angle in degree
F1y=-1e5;                       %Applied load(down at node 1)
d3x=0;                          %x-displacement 3 is zero
d3y=0;                          %y-displacement at node 3 is zero
Theta_3=0;                       %angular displacement at node 3 is zero
C1=cosd(Theta_1);                %Cosine of angle between local x and global x for element 1
S1=sind(Theta_1);                %Sine of angle between local x and global x for element 1
C2=cosd(Theta_2);                %Cosine of angle between local x and global x for element 2
S2=sind(Theta_2);                %Sine of angle between local x and global x for element 2

M1=(AA*C1^2)+((12*II*S1^2)/LL1^2)
M2=(AA-((12*II)/LL1^2))*C1*S1
M3=(AA*S1^2)+((12*II*C1^2)/LL1^2)
ISL=6*II*S1/LL1
ICL=6*II*C1/LL1
K1=(EE/LL1)*[M1 M2 -ISL -M1 -M2 -ISL;M2 M3 ICL -M2 -M3 ICL;-ISL ICL 4*II ISL -ICL 2*II;
    -M1 -M2 ISL M1 M2 ISL;-M2 -M3 -ICL M2 M3 -ICL;-ISL ICL 2*II ISL -ICL 4*II]
K1a=K1
K1a(9,9)=0


M12=(AA*C2^2)+((12*II*S2^2)/LL2^2)
M22=(AA-((12*II)/LL2^2))*(C2*S2)
M32=(AA*S2^2)+((12*II*C2^2)/LL2^2)
ISL2=6*II*S2/LL2
ICL2=(6*II*C2)/LL2
K2=(EE/LL2)*[M12 -M22 ISL2 -M12 M22 ISL2;-M22 M32 ICL2 M22 -M32 ICL2;ISL2 ICL2 4*II -ISL2 -ICL2 2*II;
    -M12 M22 -ISL2 M12 -M22 -ISL2;M22 -M32 -ICL2 -M22 M32 -ICL2;ISL2 ICL2 2*II -ISL2 -ICL2 4*II]
K2a=zeros(size(K1a))
K2a(4:end,4:end)=K2
K=K1a+K2a
Ksub=K([1,2,3,4,5,6],[1,2,3,4,5,6])
Fy=-10000
Fpart=[0;Fy;0;0;0;0]               %e partitioned vector of applied loads
[D]=linsolve(Ksub,Fpart)
syms U1 V1 Ph1 U2 V2 Ph2
U1=D(1,1)
V1=D(2,1) 
Ph1=D(3,1) 
U2=D(4,1) 
V2=D(5,1) 
Ph2=D(6,1)
U3=0
V3=0
Ph3=0
format LongE
D1=[U1;V1;Ph1;U2;V2;Ph2;U3;V3;Ph3]
FM=K*D1
AEL=AA*EE/LL1
EIL=2*EE*II/LL1
EIL1=4*EE*II/LL1
EIL2=6*EE*II/LL1^2
EIL3=12*EE*II/LL1^3
Klocal1=[AEL 0 0 -AEL 0 0;0 EIL3 EIL2 0 -EIL3 EIL2;0 EIL2 EIL1 0 -EIL2 EIL;
        -AEL 0 0 AEL 0 0;0 -EIL3 -EIL2 0 EIL3 -EIL2;0 EIL2 EIL 0 -EIL2 EIL1]
    
TCS1=[C1 S1 0 0 0 0;-S1 C1 0 0 0 0;0 0 1 0 0 0;0 0 0 C1 S1 0;0 0 0 -S1 C1 0;0 0 0 0 0 1]
F1=Klocal1*TCS1*D

%For element 2
AEL2=AA*EE/LL2
EIL2=2*EE*II/LL2
EIL12=4*EE*II/LL2
EIL22=6*EE*II/LL2^2
EIL32=12*EE*II/LL2^3
Klocal2=[AEL2 0 0 -AEL2 0 0;0 EIL32 EIL22 0 -EIL32 EIL22;0 EIL22 EIL12 0 -EIL22 EIL2;
        -AEL2 0 0 AEL2 0 0;0 -EIL32 -EIL22 0 EIL32 -EIL22;0 EIL22 EIL2 0 -EIL22 EIL12]
    
TCS2=[C2 -S2 0 0 0 0;S2 C2 0 0 0 0;0 0 1 0 0 0;0 0 0 C2 -S2 0;0 0 0 S2 C2 0;0 0 0 0 0 1]
F2=Klocal2*TCS2*[U2;V2;Ph2;U3;V3;Ph3]