function [onmax, offmax, r, pairs]=plottingcolorRFnew(onfactor, offfactor, normalmethod, back)
%normalmethod: 0, normalize respectively by the max value of each subfield; 1,
%normalized by the max value of the whole RF
%close all;
newfig=1;
enlargefactor=10;
directoryname = uigetdir('D:\li and whit labs data backup\whole cell data backedup\012609-051409\20DAQ1_28_2009');
cd(directoryname);
[xfn,xpn]=uigetfile('*.*','pick a RF file');
rfpath=fullfile(xpn,xfn);
data=load(rfpath);
sizedata=size(data);
row=sizedata(1)/2;
column=sizedata(2);
on=data(1:row,:);
off=data(row+1:2*row,:);
rfsize=size(on);
[x,y] = meshgrid(1:rfsize(2),1:rfsize(1));
%z=data;
[xi,yi] = meshgrid(1:(1/enlargefactor):rfsize(2),1:1/enlargefactor:rfsize(1));
onzoom = interp2(x,y,on,xi,yi,'bilinear');
offzoom = interp2(x,y,off,xi,yi,'bilinear');
onoffcolor=zeros([size(onzoom) 3]);
oncolor=zeros([size(onzoom) 3]);
offcolor=zeros([size(onzoom) 3]);
background=mean(mean(back));
onmax=max(max(onzoom));
offmax=max(max(offzoom));
maxval=max(onmax, offmax);
onzoom=onzoom.*(onzoom>=background);
offzoom=offzoom.*(offzoom>=background);
if normalmethod==0
    onoffcolor(:,:,1)=onzoom/onmax*onfactor;
    onoffcolor(:,:,2)=offzoom/offmax*offfactor;
    oncolor(:,:,1)=onzoom/onmax*onfactor;
    offcolor(:,:,2)=offzoom/offmax*offfactor;
else
    onoffcolor(:,:,1)=onzoom/maxval*onfactor;
    onoffcolor(:,:,2)=offzoom/maxval*offfactor;
    oncolor(:,:,1)=onzoom/maxval*onfactor;
    offcolor(:,:,2)=offzoom/maxval*offfactor;    
end
onoffcolordisp=onoffcolor;
oncolordisp=oncolor;
offcolordisp=offcolor;
Xpos=reshape(xi,prod(size(xi)),1);
Ypos=reshape(yi,prod(size(yi)),1);

if newfig
    unitwinsize=40;
    ONOFFhandle=figure;
    set(ONOFFhandle,'Position',[40 400 unitwinsize*column unitwinsize*row]);
    set(ONOFFhandle,'Name','ONOFF');
else
    figure(ONOFFhandle);
    set(ONOFFhandle,'Position',[40 400 unitwinsize*column unitwinsize*row]);
    set(ONOFFhandle,'Name','ONOFF');
end
image(Xpos,Ypos,onoffcolordisp);
axis on;

if newfig
    ONhandle=figure;
    set(ONhandle,'Position',[40 400 unitwinsize*column unitwinsize*row]);
    set(ONhandle,'Name','ON');
else
    figure(ONhandle);
    set(ONhandle,'Position',[40 400 unitwinsize*column unitwinsize*row]);
    set(ONhandle,'Name','ON');
end
image(Xpos,Ypos,oncolordisp(:,:,1)*64);
axis on;

if newfig
    OFFhandle=figure;
    set(OFFhandle,'Position',[40 400 unitwinsize*column unitwinsize*row]);
    set(OFFhandle,'Name','OFF');
else
    figure(OFFhandle);
    set(OFFhandle,'Position',[40 400 unitwinsize*column unitwinsize*row]);
    set(OFFhandle,'Name','OFF');
end
image(Xpos,Ypos,offcolordisp(:,:,2)*64);
axis on;

ButtonName=questdlg('Save figures?', 'save option', 'Yes', 'No', 'Yes');
switch ButtonName
    case 'Yes'
    dotpos=strfind(rfpath, '.');
    onsavepath=[rfpath(1:dotpos-1) 'ON' ];
    saveas(ONhandle, onsavepath, 'fig');
    offsavepath=[rfpath(1:dotpos-1) 'OFF' ];
    saveas(OFFhandle, offsavepath, 'fig');
end

onindex=find(on>background);
offindex=find(off>background);
allindex=onindex;
for i=1:length(offindex)
   if sum(offindex(i)==onindex)==0
      allindex=[allindex; offindex(i)]; 
   end
end
pairs=[on(allindex) off(allindex)];
r=corr(pairs);
r=r(1,2);
