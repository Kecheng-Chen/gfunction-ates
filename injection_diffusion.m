function [h1,h2]=injection_diffusion(value,name)
figg = figure(1);hold on;
figg2 = figure(2);hold on;
[ve1,ve2,H_array,bound,Q_array,hea_aqu,hea_cond]=deal(1,1/value(2),value(6),value(1),value(3),value(4),value(5));
syms r t0 x00 h;
H=H_array*2;
loss1=[];
% loss2=[];
h1=[];
h2=[];
year=[];
for j=1:50
    year(end+1)=bound+(j-1);
    year(end+1)=1+(j-1);
end
A=dlmread(name, '', 9, 0);
rho_c=hea_aqu*10^6;
Q=Q_array*2;
bb=0;
a=5;
b=10000000000;
lam=hea_cond;
D=lam/rho_c;
rho=rho_c;
for index=1:length(year)
ckc=A(find(A(:,2)<=H/2),:);
x0 = ckc(:,1);
y = ckc(:,2);
rangex=sort(uniquetol(x0,0.00001));
rangey=sort(uniquetol(y,0.01));
[X,Y] = meshgrid(rangex,rangey);
if rem(index,2)==0
    indexn=(index)/2;
    index2=(index)/2;
else
    indexn=50+(index+1)/2;
    index2=(index+1)/2;
end
Z = griddata(A(:,1),A(:,2),A(:,3+indexn),X,Y);
bb0=zeros(length(rangex),1);
for christ=1:length(rangex)
   B=Z(:,christ);
   bb0(christ)=trapz(rangey,B)./(H/2);
end

rangex=rangex(2:end);
bb0=bb0(2:end);
bb0_temp=rescale((bb0-14.8)./(8.8-14.8));
loca_bb0=find(bb0_temp<1e-6);
best_rangex=rangex(loca_bb0(1));
bench=trapz(rangex(1:loca_bb0(1)),bb0_temp(1:loca_bb0(1)).*rangex(1:loca_bb0(1)));

if rem(index2,5)==1
if rem(index,2)==1
   best_rangexx=best_rangex;
   color='b';
   cao=(bb0(1:loca_bb0(1))-14.8)./(8.8-14.8);
   cao=cao-min(cao);
   figure(1);
   plot(rangex(1:loca_bb0(1)),cao,color,'DisplayName',strcat('t= ',num2str(year(index)),' year'));
else
   color='m';
   cao=(bb0(1:loca_bb0(1))-14.8)./(8.8-14.8);
   cao=cao-min(cao);
   figure(2);
   plot(rangex(1:loca_bb0(1)),cao,color,'DisplayName',strcat('t= ',num2str(year(index)),' year'));
end
end

if rem(index,2)==1
h_final=pyrunfile("sol.py","i11",X=[value(1),value(2),value(3),value(4),value(5),value(6),index2]);
t=bound*60*60*24*365;
alpha=Q./(2.*pi.*H).*4.*10^6./rho_c.*ve1-D;
beta=alpha;
if index==1
alpha=Q./(4.*pi.*H).*4.*10^6./lam.*ve1;
A00=@(h,r) quadgk(@(x) exp(-x-h.*r.^2./(4.*rho_c.*D.*H)./x+(alpha-1).*log(x)-gammaln(alpha)),r.^2./(4.*D.*t),Inf,'RelTol',0,'AbsTol',1e-64);
A0=@(h,r) arrayfun(@(x) A00(h,x), r);
A1=@(h,r) arrayfun(@(x) A00(h,x), r);
A11=@(h,r) int(exp(-x00-h.*r.^2./(4.*rho_c.*D.*H)./x00+(alpha-1).*log(x00)-gammaln(alpha)),x00,r.^2./(4.*D.*t),Inf);
ref=trapz(rangex(1:loca_bb0(1)),(A0(h_final,rangex(1:loca_bb0(1)))-bb0_temp(1:loca_bb0(1))).*rangex(1:loca_bb0(1)));
else
T00=@(h,r) exp(-h.*r.^2./(2.*beta.*rho_c.*H));
T01=@(h,r) bb.*exp(1).^((-1).*h.*H.^(-1).*rho.^(-1).*t+(-1).*b.^(-1).*((sqrt( ...
  -1)*(-1)).*beta.^(1/2).*((-1).*beta.^(-1).*r.^2+2.*t).^(1/2)).^a);
T0=@(h,r) heaviside(r.^2./(2.*beta)-t).*T01(h,r)+heaviside(-r.^2./(2.*beta)+t).*T00(h,r);
T10=@(h,r) (-1).*beta.^(-2).*exp(1).^((-1/2).*beta.^(-1).*h.*H.^(-1).*r.^2.* ...
  rho.^(-1)).*h.*H.^(-2).*rho.^(-2).*t.*((-1).*h.*r.^2+beta.*H.*rho+ ...
  beta.*h.*t);
T11=@(h,r) (-1).*a.*b.^(-2).*bb.*exp(1).^((-1).*h.*H.^(-1).*rho.^(-1).*t+(-1) ...
  .*b.^(-1).*((sqrt(-1)*(-1)).*beta.^(1/2).*((-1).*beta.^(-1).*r.^2+ ...
  2.*t).^(1/2)).^a).*t.*((sqrt(-1)*(-1)).*beta.^(1/2).*((-1).* ...
  beta.^(-1).*r.^2+2.*t).^(1/2)).^a.*(r.^2+(-2).*beta.*t).^(-2).*((( ...
  -1)+a).*b.*r.^2+(-1).*a.*b.*beta.*t+(-1).*a.*((sqrt(-1)*(-1)).* ...
  beta.^(1/2).*((-1).*beta.^(-1).*r.^2+2.*t).^(1/2)).^a.*(r.^2+(-1) ...
  .*beta.*t));
T1=@(h,r) heaviside(r.^2./(2.*beta)-t).*T11(h,r)+heaviside(-r.^2./(2.*beta)+t).*T10(h,r);
T0_contp=@(h,r) 1./2.*(1-bb).*exp(-h.*t./(H.*rho_c)).*erfc((r-sqrt(2.*beta.*t))./(sqrt(D))./(2.*sqrt(t)));
T0_contn=@(h,r) -1./2.*(1-bb).*exp(-h.*t./(H.*rho_c)).*erfc(-(r-sqrt(2.*beta.*t))./(sqrt(D))./(2.*sqrt(t)));
T0_cont=@(h,r) heaviside(r.^2./(2.*beta)-t).*T0_contp(h,r)+heaviside(-r.^2./(2.*beta)+t).*T0_contn(h,r);
A0=@(h,r) max(smoothdata(T0(h,r)+T0_cont(h,r)+D.*T1(h,r),'movmedian',10),0);
A1=@(h,r) T0(h,r)+T0_cont(h,r)+D.*T1(h,r);
ref=trapz(rangex(1:loca_bb0(1)),(A0(h_final,rangex(1:loca_bb0(1)))-bb0_temp(1:loca_bb0(1))).*rangex(1:loca_bb0(1)));
end
h1(end+1)=h_final;
loss1(end+1)=abs(ref)./bench;
if rem(index2,5)==1
figure(1);
plot(rangex(1:loca_bb0(1)),A0(h1(index2),rangex(1:loca_bb0(1))),'k--');
end
else
t=(1-bound)*60*60*24*365;
alpha=Q./(2.*pi.*H).*4.*10^6./rho_c.*ve2+D;
beta=alpha;
f=@(r) A1(h_final,r);
T0=@(h,r,t0) exp(-h*t0/(H*rho_c)).*f(((r.^2+2.*alpha.*t0).^0.5));
if index==2
    f1=@(r) A11(h_final,r);
    T00=@(h,r,t0) exp(-h*t0/(H*rho_c)).*f1(((r.^2+2.*alpha.*t0).^0.5));
    temp01=matlabFunction(diff(T00,r,2));
    temp0=@(h,r,t) arrayfun(@(rr,tt) temp01(h,rr,tt),r,t);
else
    temp0=matlabFunction(diff(T0,r,2));
end
T1=@(h,r,t0) exp(h.*r.^2./(2.*beta.*H.*rho)).*...
    arrayfun(@(r1) quadgk(@(r0) -exp(-h.*r0.^2./(2.*beta.*H.*rho)).*temp0(h,r0,(r1.^2+2.*beta.*t0-r0.^2)./(2.*beta)).*r0./beta,-sqrt(r1.^2+2.*beta.*t0),r1),r);
options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',1e10);
lb = 0;
ub = [];
T0=@(h,r) max(smoothdata(max(T0(h,r,t)+D.*T1(h,r,t),0),'movmedian',20),0);
h_final=pyrunfile("sol.py","i21",X=[value(1),value(2),value(3),value(4),value(5),value(6),index2]);
% ref=trapz(rangex(1:loca_bb0(1)),T0(h_final,rangex(1:loca_bb0(1))).*rangex(1:loca_bb0(1)));
% loss2(end+1)=abs(ref-bench)./bench;
bb=T0(h_final,2.55);
pa_initial=[a,b];
flex_fun=@(pa,r) bb.*exp(-r.^(pa(1))./pa(2));
pa_final = lsqcurvefit(flex_fun,pa_initial,rangex,T0(h_final,rangex),lb,ub,options);
a=pa_final(1);
b=pa_final(2);
h2(end+1)=h_final;
if rem(index2,5)==1
figure(2);
plot(rangex(1:end),T0(h_final,rangex(1:end)),'g--');
hold on;
end
end
end
figure(1);
xlim([0,best_rangexx]);
ylim([0,1]);
xlabel('r');
ylabel('T^*');
legend('Benchmark','Analytical');
title('Dimensionless temperature profile after injection');
figure(2);
xlim([0,best_rangexx]);
ylim([0,1]);
xlabel('r');
ylabel('T^*');
legend('Benchmark','Analytical');
title('Dimensionless temperature profile after extraction');
saveas(figg,strcat("results/profile_after_injection.png"));
saveas(figg2,strcat("results/profile_after_extraction.png"));
disp(strcat('Mean energy estimation error rate after each injection: ',num2str(mean(loss1)*100),'%'));
end