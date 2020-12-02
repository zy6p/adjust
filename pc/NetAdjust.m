function  NetAdjust(T1,T21,T22,T31,T32,T33,sigma1,sigma2,num1,num2)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
for ii=1:9
    if T1(ii,1)==0
        T1(ii,1)=9;
    end
    if T1(ii,2)==0
        T1(ii,2)=9;
    end
end
T11=T1(:,1)*10+T1(:,2);
T12=T1(:,3);
global cX cY cD
L0=[T22;T12];
X0=zeros(2*(num1-length(T31)),1);
cX=zeros(num1,1);
cY=zeros(num1,1);
cD=1:num1;
cD(T31)=0;
cX(T31)=T32;
cY(T31)=T33;
ox=[269.011;163.475;152.66;104.35;172.45;278.44;164.65;193.77;284.43];
oy=[292.728;292.866;264.52;183.83;123.89;124.98;253.76;209.33;242.22];
for ii=1:num1
    if cD(ii)>0
        if ii>T31(1)
            cD(ii)=cD(ii)-1;
        end
        if ii>T31(2)
            cD(ii)=cD(ii)-1;
        end
    end
end
for ii=1:num1
    if cD(ii)~=0    
        cX(ii)=ox(ii);
        cY(ii)=oy(ii);
        X0(2*cD(ii)-1)=cX(ii);
        X0(2*cD(ii))=cY(ii);
    end
end
P=[eye(length(T21)),zeros(length(T21),length(T11));zeros(length(T11),length(T21)),(sigma1/206265)^2*eye(length(T11))/(diag(T12)*sigma2)^2];
x=ones(2*(num1-length(T31)),1);
B=zeros(length(T21)+length(T11),2*(num1-length(T31)));
l=zeros(length(T21)+length(T11),1);
iabs=0;
while max(abs(x))>1e-100
for ii=1:length(T21)
    B(ii,:)=ffoc(T21(ii));
    l(ii)=L0(ii)-ttoc(T21(ii));
end
for ii=length(T21)+1:length(T21)+length(T11)
    B(ii,:)=ssoc(T11(ii-length(T21)));
    l(ii)=L0(ii)-rroc(T11(ii-length(T21)));
end
    x=(B'*P*B)\B'*P*l;
    X0=X0+x;
    for jj=1:num1
        if cD(jj)~=0
            cX(jj)=X0(2*cD(jj)-1);
            cY(jj)=X0(2*cD(jj));
        end
    end
    iabs=iabs+1
end
v=B*x-l;
Lg=L0+v;
vg=[3600*rad2deg(v(1:length(T21)));1000*v(length(T21)+1:end)];
disp('边观测值平差结果\n边编号\t观测值（m）\t平差值\t改正数（m）')
vpa([T11 T12 Lg(length(T21)+1:end) v(length(T21)+1:end)],6)
Lq=rad2deg(Lg(1:length(T21)));
Lw=fix(Lq/1);
Le=fix(60*(Lq-Lw));
Lr=(Lq-Lw)*3600-60*Le;
Lq1=rad2deg(L0(1:length(T21)));
Lw1=fix(Lq1/1);
Le1=fix(60*(Lq1-Lw1));
Lr1=(Lq1-Lw1)*3600-60*Le1;
disp('角度观测值平差结果\n编号\t角观测值（°\t′\t″）\t角平差值（°\t′\t″）\t角改正数（″）')
vpa([T21 Lw1 Le1 Lr1 Lw Le Lr vg(1:length(T21))],6)
sigma5=sqrt(v'*P*v/8);
% sigma6=sqrt(vg'*P*vg/8);
DXX=eye(2*(num1-length(T31)))/(B'*P*B)*sigma5*sigma5;
Dsigma=zeros(num1,1);
QLL=B/(B'*P*B)*B';
QSS=QLL(length(T21)+1:length(T21)+length(T11),length(T21)+1:length(T21)+length(T11));
sigma3=sigma5*sqrt(max(max(QSS)));
for ii=1:num1
    if cD(ii)~=0
        Dsigma(ii)=1000*sqrt(DXX(2*cD(ii)-1,2*cD(ii)-1)+DXX(2*cD(ii),2*cD(ii)));
    end
end
sigma4=max(Dsigma);
dx=zeros(num1,2);
    for jj=1:num1
        if cD(jj)~=0
            dx(jj,1)=x(2*cD(jj)-1);
            dx(jj,2)=x(2*cD(jj));
        end
    end
fprintf('平差后最弱边相对中误差%f\n平差后最弱点位中误差%fmm\n平差后单位权中误差=%f″',sigma3,sigma4,206265*sigma5)
disp('坐标值平差结果\n编号\tX平差值（m）\tY平差值（m）\tX改正数（m）\tY改正数（m）\t点位中误差（mm）')
vpa([cX,cY,dx,Dsigma],6)

for jj=1:num1
    uuoc(T11(jj));
end
plot(cY(T31),cX(T31))
axis equal

    function uuoc(aa)
        cc1=fix(aa/10);
        cc2=aa-cc1*10;
        plot([cY(cc1),cY(cc2)],[cX(cc1),cX(cc2)])
        hold on
    end

    function [ddd]=ttoc(aaa)
        cc1=fix(aaa/100);
        cc2=fix((aaa-cc1*100)/10);
        cc3=aaa-cc1*100-cc2*10;
        ddd=angle(cX(cc1),cY(cc1),cX(cc2),cY(cc2),cX(cc3),cY(cc3));
        function [oangle] = angle(a,b,c,d,e,f)
            function [ccan]=cccw(p,q,m,n)
                if n-q>0 && m-p>0
                    ccan=atan((n-q)/(m-p));
                elseif n-q>0 && m-p<0
                    ccan=pi+atan((n-q)/(m-p));
                elseif n-q<0 && m-p>0
                    ccan=2*pi+atan((n-q)/(m-p));
                else
                    ccan=pi+atan((n-q)/(m-p));
                end
            end
            oangle=cccw(c,d,e,f)-cccw(a,b,c,d)+pi;
            if oangle<0
                oangle=oangle+2*pi;
            elseif oangle>2*pi
                oangle=oangle-2*pi;
            end
        end
    end

    function [dd]=rroc(aa)
        cc1=fix(aa/10);
        cc2=aa-cc1*10;
        dd=sqrt((cX(cc1)-cX(cc2))^2+(cY(cc1)-cY(cc2))^2);
    end


    function [ddd]=ffoc(aaa)
        cc1=fix(aaa/100);
        cc2=fix((aaa-cc1*100)/10);
        cc3=aaa-cc1*100-cc2*10;
        bb1=cD(cc1);
        bb2=cD(cc2);
        bb3=cD(cc3);
        ddd=zeros(1,2*(num1-length(T31)));
        if bb1~=0
            ddd(2*bb1-1)=(cY(cc1) - cY(cc2))/(((cY(cc1) - cY(cc2))^2/(cX(cc1) - cX(cc2))^2 + 1)*(cX(cc1) - cX(cc2))^2);
            ddd(2*bb1)= -1/(((cY(cc1) - cY(cc2))^2/(cX(cc1) - cX(cc2))^2 + 1)*(cX(cc1) - cX(cc2)));
        end
        if bb2~=0
            ddd(2*bb2-1)=- (cY(cc1) - cY(cc2))/(((cY(cc1) - cY(cc2))^2/(cX(cc1) - cX(cc2))^2 + 1)*(cX(cc1) - cX(cc2))^2) - (cY(cc2) - cY(cc3))/(((cY(cc2) - cY(cc3))^2/(cX(cc2) - cX(cc3))^2 + 1)*(cX(cc2) - cX(cc3))^2);
            ddd(2*bb2)=1/(((cY(cc1) - cY(cc2))^2/(cX(cc1) - cX(cc2))^2 + 1)*(cX(cc1) - cX(cc2))) + 1/(((cY(cc2) - cY(cc3))^2/(cX(cc2) - cX(cc3))^2 + 1)*(cX(cc2) - cX(cc3)));
        end
        if bb3~=0
            ddd(2*bb3-1)= (cY(cc2) - cY(cc3))/(((cY(cc2) - cY(cc3))^2/(cX(cc2) - cX(cc3))^2 + 1)*(cX(cc2) - cX(cc3))^2);
            ddd(2*bb3)=-1/(((cY(cc2) - cY(cc3))^2/(cX(cc2) - cX(cc3))^2 + 1)*(cX(cc2) - cX(cc3)));
        end
    end


    function [dd] = ssoc(aa)
        %UNTITLED6 此处显示有关此函数的摘要
        %   此处显示详细说明
        cc1=fix(aa/10);
        cc2=aa-cc1*10;
        bb1=cD(cc1);
        bb2=cD(cc2);
        dd=zeros(1,2*(num1-length(T31)));
        if bb1~=0
            dd(2*bb1-1)=(2*cX(cc1) - 2*cX(cc2))/(2*((cX(cc1) - cX(cc2))^2 + (cX(cc1) - cX(cc2))^2)^(1/2));
            dd(2*bb1)=(2*cY(cc1) - 2*cY(cc2))/(2*((cX(cc1) - cX(cc2))^2 + (cX(cc1) - cX(cc2))^2)^(1/2));
        end
        if bb2~=0
            dd(2*bb2-1)=(2*cX(cc2) - 2*cX(cc1))/(2*((cX(cc1) - cX(cc2))^2 + (cX(cc1) - cX(cc2))^2)^(1/2));
            dd(2*bb2)=(2*cY(cc2) - 2*cY(cc1))/(2*((cX(cc1) - cX(cc2))^2 + (cX(cc1) - cX(cc2))^2)^(1/2));
        end
    end
end