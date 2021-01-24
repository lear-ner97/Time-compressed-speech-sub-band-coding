clear all
clc
close all
%%seance 1
%1**********
[s,fe]=audioread('voix_homme_8.wav');
soundsc(s);
Lfft=1024;
%2**********
p=24;fc=fe/4;
fcn=fc/(fe/2);
h=fir1(p,fcn,rectwin(p+1)); %filtre passe-bas d'ordre p et de freq de 
%coupure fc / rectwin(x)==> x quelconque

%3gain amplitude**************
[H,f]=freqz(h,1,1000,fe);  %freqz(a,b,fe)   
                        % H est la reponse frequentielle et f est le 
                        %vecteur freq
 figure(1);plot(f,abs(H));hold on;grid on;
 %commentaire : présente des ondulations importantes

 %4:choisir une fenetre de hamming***********
hh=fir1(p,fcn,hamming(p+1));
[HB,f]=freqz(hh,1,1000,fe);
figure(1);plot(f,abs(HB));legend('fenêtre rectangulaire','fenêtre de Hamming');grid on;
%contrairement à la fenetre rectangulaire, la fenetre de hamming
%ne présente pas d'ondulations 

%5*************
% pour p=24 la bande de transition est plus etroite ===> on augmente p 
% pour réduire la bande de transition


%spectre de notre signal audio(question supplem)
Lfft=1024;
N=length(s);
f1= (-Lfft/2 : Lfft/2-1)*fe/Lfft;
xs=fftshift(fft(s,Lfft));figure(2);subplot(2,1,1);
plot(f1,abs(xs)); 
grid on; 
xlabel('Frequence (Hz)'); 
ylabel('Amplitude'); title('Spectre d’amplitude du signal s') ;

%choix : *on augmente l'ordre du filtre p pour réduire la bande de transition 
         %(wp diminue)
         %ET
       % *on choisit la fenetre de hamming pour réduire les ondulations

%6************
s1=filter(h,1,s);%filter(ha,b,x) appliquée dans le domaine temporel : s1 
%temporel et s temporel
S1=fftshift(fft(s1,Lfft));
figure(2);subplot(2,1,2);plot(f1,abs(S1));
grid on;title('spectre du signal filtré par le passe-bas');

%%seance2

%7 realiser un filtre passe-haut***************
high=fir1(p,fcn,'high',hamming(p+1)); % c la bande Haute Frequence du signal
%hamming(p+1) est facultative (g imposé cette fenetre sinon il choisit une
%fenetre de son choix)
[HH,f]=freqz(high,1,1000,fe);
figure(3);plot(f,abs(HH));title('Réponse fréquentielle du filtre passe-haut');grid on;
%sa reponse frequencielle
HS=HH+HB;figure(4);plot(f,abs(HS));title('somme filtre H+B (Hamming)');grid on;
%somme filtre haut+filtre bas(hamming)

%8***************
%représentation 
%N=length(s);
f1= (-Lfft/2 : Lfft/2-1)*fe/Lfft;
xs=fftshift(fft(s,Lfft));figure(5);subplot(2,1,1);plot(f1,abs(xs));
title('spectre de s');grid on;

s2=filter(high,1,s);
S2=fftshift(fft(s2,Lfft));
figure(5);subplot(2,1,2);plot(f1,abs(S2));grid on;
title('spectre de s filtré par un filtre passe-haut');

%9***************
sr=s1+s2; % = s (on récupère le signal s)
figure(6);subplot(2,2,1);plot(s);title('signal s');grid on;
subplot(2,2,2);plot(sr);title('somme s1+s2');grid on;
soundsc(sr); 

%%Seance3:
%*****************************DECIMATION*****************************
%1*******
M=2;fe1=fe/2;
s1d=s1([1:M:length(s1)]);
s2d=s2([1:M:length(s2)]);


%2********

freq1=(-Lfft/2:Lfft/2-1)*fe1/Lfft;
Yf1=fftshift(fft(s1d,Lfft));
Yf2=fftshift(fft(s2d,Lfft));

figure(7);subplot(2,2,1);plot(freq1,abs(Yf1));title('s1d');grid on;
subplot(2,2,2);plot(freq1,abs(Yf2));title('s2d'); grid on;

%//////////////////////////Quantification/////////////////////////

%1****************
%function [xq,rsb]=numerise(x,l);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% La fonction numerise effectue la quantification uniforme d'un signal "x" 
% en utilisant "l" bits et calcule le signal résultant "xq" ainsi que le 
% rapport signal sur bruit "rsb" correspondant.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%L=2^l;                         % Nbre de niveaux de quantifications 
%M = max(abs(x));

%xq = (2*M/L).*round(L.*x./(2*M));
%rsb = 20*log10(std(x)/std(x-xq));
%Appliquons cela sur les deux signaux s1 et s2

%2************************
l=2:16;
rsb1=zeros(1,15);rsb2=zeros(1,15);
for i=2:16
   [sq1,rb1]=numerise(s1d,i);
   rsb1(i-1)=rb1;
   [sq2,rb2]=numerise(s2d,i);
   rsb2(i-1)=rb2;
 
end
figure(8);plot(l,rsb1);hold on;plot(l,rsb2);hold off;
title('Quantification');legend('sd1','sd2');grid on;

%je choisis l1=8 et l2=12 sur lesquels je vais travailler par la suite

[sq1,rsb1]=numerise(s1d,7);
[sq2,rsb2]=numerise(s2d,13);
%Interpolation\\\\\\\\\\\\\\\\\\\\\\\\\\

%1*****************
M=2;
L1=length(sq1);
y1=zeros(1,M*L1);         y2=zeros(1,M*L1);
y1([1:M:length(y1)])=sq1;         y2([1:M:length(y2)])=sq2;


freq1=(-Lfft/2:Lfft/2-1)*fe/Lfft; %voir doc: et ce en ramenant la fréq d'éch à Fe
Yf1=fftshift(fft(y1,Lfft));
Yf2=fftshift(fft(y2,Lfft));

%2***********
figure(9);subplot(2,2,1);plot(freq1,abs(Yf1));
title('sq1-suréchantillonnage');grid on;
subplot(2,2,2);plot(freq1,abs(Yf2));title('sq2-suréchantillonnage');
grid on;

%Filtres de synthèse////////////////

%1***********
ss1=filter(h,1,y1);
SS1=fftshift(fft(ss1,Lfft));
figure(10);subplot(2,1,2);plot(f1,abs(SS1));grid on;
title('y1-filtrage BF');%voir title

ss2=filter(high,1,y2);
SS2=fftshift(fft(ss2,Lfft));
figure(10);subplot(2,2,2);plot(f1,abs(SS2));grid on;
title('y2-filtrage HF');%voir title

%2************
figure(11);plot(f1,abs(SS2+SS1));%spectre de la somme des 2 signaux sur-echan et filtrés
title('somme des spectres de y1 et y2');grid on;
soundsc(ss1+ss2);%sounds like s


%décalage temporel résultant=(l1+l2)*p

%3
%On varie les valeurs de l afin d'obtenir le
%signal le plus proche du signal de départ.

%4
%Débit binaire résultant de toute la chaine selon les l choisi
%D1=l1*fe;
%D2=l2*fe;
%Dr=(D1+D2)/2

p=length(ss1+ss2);

figure(12),plot(ss1+ss2)

hold on
plot([0:N-1],s,'r'),legend('Signal interpolé','Signal original');

%D=Nbr_sig*nbr_bit*fe
%%seance4
