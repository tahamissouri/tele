Fe=24000;
Rb=3000;
n=8000;
%Modulateur 1 :
Ts1=1/Rb;
nb=Ts1*Fe;
bits=randi([0,1],1,n);
bits1=2*bits-1;

nrzb=zeros(n,nb);
nrzb(:,1)=bits1;
nrzb=reshape(nrzb.',1,[]);
%plot(0:1/Fe:n*nb/Fe-1/Fe,nrzb);
h1=ones(1,nb);
%{
fnrzb=filter(h1,1,nrzb);
figure(1)
plot(fnrzb);
DSP1=pwelch(fnrzb,[],[],Fe,'twosided');
figure(2)
hold on
semilogy(fftshift(DSP1));

f=linspace(-Fe/2,Fe/2,length(DSP1));
DSP1th=0.25*Ts1*sinc(f*Ts1).^2;
semilogy((DSP1th));
hold off
%}
%%Modulateur2:
M=4;
RS2=Rb/log2(M);
NS2=Fe/RS2;
H2=ones(1,NS2);

bits2=randi([0,1],1,n);
sym=(2*bi2de(reshape(bits2,2,length(bits2)/2).')-3);
dirac=kron(sym,[1 zeros(1,NS2-1)]);
signal=filter(H2,1,dirac);
DSP2=pwelch(signal,[],[],[],Fe,'twosided');

figure(4)
hold on
semilogy(fftshift(DSP2));

f=linspace(-Fe/2,Fe/2,length(DSP2));
DSP2th=0.25*Ts1*sinc(f*Ts1).^2;
semilogy((DSP2th));
hold off

%% modulateur3:
H3=rcosdesign(0.5,8,Fe/Rb-1);
kr=kron(2*bits2-1,[1 zeros(1,Fe/Rb)]);
fsignal3=filter(H3,1,kr);

DSP3=pwelch(fsignal3,[],[],[],Fe,'twosided');
figure(6)
hold on
semilogy(fftshift(DSP3));
y=0:1/Fe:n/Fe;
%semilogy(fftshift(fx(y)));
%function DfnrzbSP3th=fx(x)
 % if abs(x)<=0.5/2*Ts1
  %    DSP3th=1;
  %elseif  0.5/2*Ts1<=abs(x) && abs(x)<=1.5/2*Ts1
   % DSP3th=0.5*(1+cos(0.5*pi*Ts1(abs(x)-0.5/2*Ts1)));
  %else
   % DSP3th=0;
  %end
%end


%% Etude des interf ́erences entre symbole et du critere de Nyquist
NS1=8;
g=conv(h1,h1);
signal_sortie=filter(g,1,nrzb);
figure(1)
plot(g);
%figure(2)
%plot(signal_sortie);
%figure(3)
%plot(reshape(signal_sortie,NS1,length(signal_sortie)/NS1))
%r=TEB(NS1,signal_sortie,NS1,bits);
%%
BW=8000;
nr=61;
hc=2*(BW/Fe)*sinc(2*(BW/Fe)*(-nr/2:1/Fe:nr/2));
chaine_trans=conv(g,hc);
signal_recu=filter(hc,1,signal_sortie);
figure(3)
plot(chaine_trans);
figure(4)
plot(reshape(signal_recu,NS1,length(signal_recu)/NS1))
figure(5)
semilogy(abs(fft(g)));
hold on
semilogy(abs(fft(hc)));
hold off
figure(7)
plot(signal_recu);
%% partie 4
% 1-
sigma=0.22;

bruit=sigma*randn(1,length(nrzb));
x_1=filter(h1,1,nrzb) ;
y_1=filter(hc,1,x_1) ;
h2=conv(h1,hc);
h3=conv(H2,hc);
signal_bruitee1 =y_1+bruit ;
x_2=filter(H2,1,nrzb) ;
y_2=filter(hc,1,x_2) ;
signal_bruitee2 =y_2+bruit ;

signal_sortie1=filter(h1,1,signal_bruitee1);
signal_sortie2=filter(h1,1,signal_bruitee2);

figure(8)
plot(signal_sortie1);
hold on
plot(signal_sortie2);
hold off

TEB1=TEB(NS1,signal_sortie1,NS1,bits) ;
TEB2=TEB(NS1,signal_sortie2,NS1,bits) ;
t=NS1:NS1:NS1*n ;
signal_ech1 = signal_sortie1(t) ;
signal_ech2 = signal_sortie2(t);

figure(10)
plot(signal_ech1);
hold on
plot(signal_ech2);
hold off

snr_db =0:1:10 ;
bTEB =zeros (length(snr_db),1);
for i=1:1:10
    snr = 10^(i/10);
    bbits=randi([0,1],1,n);
    bbits1=2*bbits-1;
    bnrzb=zeros(n,nb);
    bnrzb(:,1)=bbits1;
    bnrzb=reshape(bnrzb.',1,[]);
    bfnrzb=filter(h1,1,bnrzb);
    n0=((mean(abs(bfnrzb).^2))*NS1)/(2*snr);
    b=sqrt(n0)*randn(1,length(bnrzb));
    r=bfnrzb+b;
    b2=filter(hc,1,r);
    bTEB(i) = TEB(NS1,b2,NS1,bbits);
end
figure(11)
semilogy(snr_db,bTEB);
% 2-















