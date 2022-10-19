function [h1,h2]=injection_advection(value,name)
figg = figure(1);hold on;
figg2 = figure(2);hold on;
[ve1,ve2,H_array,bound,Q_array,hea_aqu,hea_cond]=deal(1,1/value(2),value(6),value(1),value(3),value(4),value(5));
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
x=0.1:1:2000;
rho_c=hea_aqu*10^6;
Q=Q_array*2;
bb=0;
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
h_final=pyrunfile("sol.py","i10",X=[value(1),value(2),value(3),value(4),value(5),value(6),index2]);
t=bound*60*60*24*365;
alpha=Q./(2.*pi.*H).*4.*10^6./rho_c.*ve1-D;
beta=alpha;
if index==1
alpha=Q./(4.*pi.*H).*4.*10^6./lam.*ve1;
A00=@(h,r) quadgk(@(x) exp(-x-h.*r.^2./(4.*rho_c.*D.*H)./x+(alpha-1).*log(x)-gammaln(alpha)),r.^2./(4.*D.*t),Inf,'RelTol',0,'AbsTol',1e-64);
A0=@(h,r) arrayfun(@(x) A00(h,x), r);
ref=trapz(rangex(1:loca_bb0(1)),(A0(h_final,rangex(1:loca_bb0(1)))-bb0_temp(1:loca_bb0(1))).*rangex(1:loca_bb0(1)));
else
% f=@(r) T0(h_final,r);
f=@(r) arrayfun(@(r0) feval(SplineFit,r0),r);
T01=@(h,r) exp(1).^((-1).*h.*H.^(-1).*rho.^(-1).*t).*f(real((sqrt(-1)*(-1)).* ...
  beta.^(1/2).*((-1).*beta.^(-1).*r.^2+2.*t).^(1/2)));
T00=@(h,r) exp(-h.*r.^2./(2.*beta.*rho_c.*H));
T0_contp=@(h,r) 1./2.*(1-bb).*exp(-h.*t./(H.*rho_c)).*erfc((r-sqrt(2.*beta.*t))./(sqrt(D))./(2.*sqrt(t)));
T0_contn=@(h,r) -1./2.*(1-bb).*exp(-h.*t./(H.*rho_c)).*erfc(-(r-sqrt(2.*beta.*t))./(sqrt(D))./(2.*sqrt(t)));
T0_cont=@(h,r) heaviside(r.^2./(2.*beta)-t).*T0_contp(h,r)+heaviside(-r.^2./(2.*beta)+t).*T0_contn(h,r);
A0=@(h,r) heaviside(r.^2./(2.*beta)-t).*T01(h,r)+heaviside(-r.^2./(2.*beta)+t).*T00(h,r)+T0_cont(h,r);
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
seq=0.1:0.1:2000;
f=@(r) A0(h_final,r);
best_ckc=f(seq);
loca=find(best_ckc<1e-10);
best_r=seq(loca(1));
fun10=@(h,r) exp(-h.*t./(H.*rho)).*f(sqrt(r.^2+2.*beta*t));
T0=@(h,r) real(arrayfun(@(r0) fun10(h,r0),r));
h_final=pyrunfile("sol.py","i20",X=[value(1),value(2),value(3),value(4),value(5),value(6),index2]);
SplineFit = fit(x(1:end)',max(T0(h_final,x(1:end))',0), 'smoothingspline');
% ref=trapz(x(1:end)',feval(SplineFit,x(1:end)').*x(1:end)');
% loss2(end+1)=abs(ref-bench)./bench;
h2(end+1)=h_final;
bb=feval(SplineFit,2.55);
if rem(index2,5)==1
figure(2);
plot(x(1:end)',feval(SplineFit,x(1:end)'),'g--');
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