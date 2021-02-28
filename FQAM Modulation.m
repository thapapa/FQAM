clc;
clear all;
close all;


M_Q=input(' enter the M_Q value : ');
M_F=input(' enter the M_F value : ');
Oversamping_Factor=input(' enter oversampling factor : ');

fprintf('\n\n\n');

Eb=1;                 %Energy Per Bit                     
SNRdB=-5:1:10;        %SNR in dB
SNR=10.^(SNRdB/10);   
ber =[];

M= M_F*M_Q;           %Modulation Order
Ld=log2(M);
ds=ceil(Ld);
dif=ds-Ld;
bp=.000001;           %bit period
k=1;
if(dif~=0)
   error('the value of M is only acceptable if log2(M)is an integer');
end
%XXXXXXXXXXXXXXXXXXX binary Information Generation XXXXXXXXXXXXXXXXXXXXXXXX
nbit=input("Enter number of bit nbit:");
nbit=nbit-(mod(nbit,log2(M)));    %nbit proper for matrix ,(nbit/log2(M) x log2(M))

                                      
msg=round(rand(nbit,1));          % information generation as binary form
disp(' binary information at transmitter ');
disp(msg);
fprintf('\n\n');
 


% binary information convert into symbolic form for M-array FQAM modulation

M=M;                  % order of FQAM modulation
msg_reshape=reshape(msg,log2(M),nbit/log2(M))'; %Reshape to Symbols
disp(' information are reshaped for convert symbolic form');
disp(msg_reshape);
fprintf('\n\n');
size(msg_reshape);



for(j=1:1:nbit/log2(M))
   for(i=1:1:log2(M_F))
       a(j,i)=num2str(msg_reshape(j,i));  %Frequency Bits    
   end
   
   for(i=(log2(M_F)+1):1:log2(M))
       b(j,k)=num2str(msg_reshape(j,i));   %QAM bits
       k=k+1;
   end
  k=1; 
end

Qas=bin2dec(b);         %QAM Bits Binary to Decimal
Qass=Qas';
Fas=bin2dec(a);         %Frequency Bits Binary to Decimal
Fass=Fas';

p=qammod(Qass,M_Q);     %constalation design for M-array QAM acording to symbol
sym=0:1:M_Q-1;     % considerable symbol of M-array QAM, just for scatterplot
pp=qammod(sym,M_Q);                     %constalation diagram for M-array QAM         
scatterplot(pp),grid on;
title('consttelation diagram for M-array QAM');




%XXXXXXXXXXXXXXXXXXXXXX  FQAM modulation XXXXXXXXXXXXXXXXXXXXXXXXXXX
             %%%%%%%  Transmitter   %%%%%%%  
             
             
RR=real(p);                %QAM Real Part 
II=imag(p);                %QAM Imaginary Part
sp=bp*log2(M);             %Symbol Period
sr=1/sp;                   % symbol rate
f_start=sr;                %Starting Frequency 

t=sp/Oversamping_Factor:sp/Oversamping_Factor:sp;    %n=T/OSF
ss=length(t);
m=[];


                  %Normalizing factor 
if M_Q == 4
 a= 1/sqrt(2);
else if M_Q == 16                           
        a = sqrt(1/10);
     end
 
end

                %Transmission      
 for(i=1:1:length(Fass))

  f=f_start*(Fass(i)+1);           %Frequency Selection
   
  yr=a*RR(i)*cos(2*pi*f*t);                     %  real component
  yim=a*II(i)*sin(2*pi*f*t);            % Quadrature  component 
  y=yr+yim;
  
 
  m=[m y];
 end

                 %Plot FQAM Signal
tt=sp/Oversamping_Factor:sp/Oversamping_Factor:sp*length(Fass);
figure(1);
subplot(3,1,3);
plot(tt,m);
title('waveform for MFQAM modulation' );
xlabel('time(sec)');
ylabel('amplitude(volt)');




     %%%%%%%%%%%%%%%%%%%%%   AWGN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%   Frequency Trap    %%%%%%%%%
         %In trap table there is the central frequency id for every 
         %M_F Frequencies
trap=[];
for(i=1:1:M_F)
f=f_start*i
   
  kr=cos(2*pi*f*t);                     
  kim=sin(2*pi*f*t);             
  k=kr+kim;
  h = abs(fft(k));
  h(1) = [];
  tmax = max(h);
  h(ceil(end/2):end) = [];
  abovecutoff = h > tmax / 2;   %3 dB is factor of 2
  lowbin  = find(abovecutoff, 1, 'first');
  highbin = sum(abovecutoff);
  centbin = sqrt(lowbin * highbin);
  trap=[trap centbin];
end


               %%%%%%%%%%%%  Receiver   %%%%%%%%%%
for count=1:1:length(SNR)    
      snumber=0;
      fnumber=0;
      qnumber=0;
      bother=0;
      tro=[];
      

No=Eb/SNR(count);    
Noise=sqrt(No/2)*(randn(1,(ss*length(Fass)))+1i*randn(1,(ss*length(Fass)))); 
m_new= m + Noise;   %Noise added Signal


if count==5
    tt=sp/Oversamping_Factor:sp/Oversamping_Factor:sp*length(Fass);
    figure(2);
    subplot(3,1,3);
    plot(tt,m_new);
    title('waveform for MFQAM modulation');
    xlabel('time(sec)');
    ylabel('amplitude(volt)');
end


                
m_reshape=reshape(m_new,Oversamping_Factor,(length(m_new)/Oversamping_Factor))';  %%symbol reshape
l=(length(m_new)/Oversamping_Factor);
V=[];
clear j;

for(i=1:1:l)
          %%%   Frequency Carrier Detection
h = abs(fft(m_reshape(i,:)));
h(1) = [];
tmax = max(h);
h(ceil(end/2):end) = [];
abovecutoff = h > tmax /3;   %3 dB is factor of 2
lowbin  = find(abovecutoff, 1, 'first');
highbin = sum(abovecutoff);
centbin = sqrt(lowbin * highbin);

    A(1:M_F)=-centbin;
    
    B=trap+A; 
  
    minimum=min(min(abs(B)));
    [x]=find(abs(B)==minimum);           %%% x=active carrier
    
    q=x-1;                               %% M-F ary bits
    tro=[tro q];
    p=de2bi(q,log2(M_F),'left-msb');
    
    f=f_start*x;
  %%%% QAM symbols 
      y1=cos(2*pi*f*t);                                     % inphase component
      y2=sin(2*pi*f*t);                                  % quadrature component
      mm1=y1.*(m_reshape(i,:)./a); 
      mm2=y2.*(m_reshape(i,:)./a);                                    
      z1=trapz(t,mm1);                                             % integration
      z2=trapz(t,mm2);                                            % integration
      zz1=round(2*z1/sp);
      zz2=round(2*z2/sp);
      gt=zz1+j*zz2;
      ax=qamdemod(gt,M_Q);
      bi_in=de2bi(ax,log2(M_Q),'left-msb');          %%% M_Q ary bits
      V=[V p bi_in];                               %%%% Reconstruction
end
tf = isequal(msg',V);                          
[number,ratio] = biterr(msg',V);


ber(count)=ratio;
diffelements=sum(Fass~=tro);
result=diffelements/numel(Fass);
%%%%%%%%%%%%SER

V_reshape=reshape(V,log2(M),nbit/log2(M))';
  for i=1:1:l
    sf=isequal(msg_reshape(i,:),V_reshape(i,:));

    if sf==0
      snumber= snumber + 1;
      ff=isequal(msg_reshape(i,(1:log2(M_F))),V_reshape(i,(1:log2(M_F))));
      qf=isequal(msg_reshape(i,(log2(M_F)+1):log2(M)),V_reshape(i,(log2(M_F)+1):log2(M)));
      if ff==0%(msg_reshape(i,(1:log2(M_F))) ~= V_reshape(i,(1:log2(M_F))))
          fnumber=fnumber+1;
      end
      if qf==0%(msg_reshape(i,(log2(M_F)+1):log2(M)) ~= V_reshape(i,(log2(M_F)+1):log2(M)))
          qnumber=qnumber+1;
      end
      if ff==0 && qf==0
          bother = bother +1;
      end
      %if ((msg_reshape(i,(1:log2(M_F)))~=V_reshape(i,(1:log2(M_F))))&&(msg_reshape(i,(log2(M_F)+1):log2(M))~=V_reshape(i,(log2(M_F)+1):log2(M))))==1
       %   dnumber ++;
      %end
    end
  end

ser(count)=snumber/l;
fer(count)=fnumber/snumber;
qer(count)=qnumber/snumber;
aller(count)=bother/snumber;
%der(count)=dnumber/snumber;
%%%%%frequency error count
fnumber
qnumber
bother

%%%%%qam error count

end
berQ4 = berawgn(SNRdB,'qam',4);
berQ16 = berawgn(SNRdB,'qam',16);
berQ64 = berawgn(SNRdB,'qam',64);
 


figure(3);

semilogy(SNRdB,ber,'-bo',SNRdB,berQ4,'-mx',SNRdB,berQ16,'-*r',SNRdB,berQ64,'-.k')
hold on;
grid on;
 legend('Fqam M_Q= 4 M_F = 4','4 qam','16 qam','64 qam');
 axis([min(SNRdB) max(SNRdB) 10^(-5) 1]);
 title('BER Vs SNR FQAM OVERSAMLPING 100');
 xlabel('SNRdB'); ylabel('BER');
 hold off
























