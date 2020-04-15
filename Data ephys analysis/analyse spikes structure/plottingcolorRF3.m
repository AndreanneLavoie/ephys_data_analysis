function maxval=plottingcolorRF3(onfactor, offfactor, normalmethod, back)

% for plotting of 1x15 RF
close all;
newfig=1;
enlargefactor=10;
directoryname = uigetdir('E:\extarcellular\Andreanne\code\');
cd(directoryname);
[xfn,xpn]=uigetfile('*.*','pick a RF file');
rfpath=fullfile(xpn,xfn);
data=load(rfpath);
sizedata=size(data);
row=sizedata(1)/2;
column=sizedata(2);
on=data(1:row,:);
off=data(row+1:2*row,:);
on=[on;on];
off=[off; off];
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
onmax=max(max(onzoom-background));
offmax=max(max(offzoom-background));
maxval=max(onmax, offmax);
if normalmethod==0
    onoffcolor(:,:,1)=(onzoom-background)/onmax*onfactor;
    onoffcolor(:,:,2)=(offzoom-background)/offmax*offfactor;
else
    onoffcolor(:,:,1)=(onzoom-background)/maxval*onfactor;
    onoffcolor(:,:,2)=(offzoom-background)/maxval*offfactor;
end
onoffcolordisp=onoffcolor.*(onoffcolor>=0);
if normalmethod==0
    oncolor(:,:,1)=(onzoom-background)/onmax*onfactor;
    offcolor(:,:,2)=(offzoom-background)/offmax*offfactor;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    oncolor(:,:,1)=(onzoom-background)/maxval*onfactor;
    offcolor(:,:,2)=(offzoom-background)/maxval*offfactor;
end
oncolordisp=oncolor.*(oncolor>=0);
offcolordisp=offcolor.*(offcolor>=0);
Xpos=reshape(xi,prod(size(xi)),1);
Ypos=reshape(yi,prod(size(yi)),1);

if newfig
    unitwinsize=40;
    ONOFFhandle=figure;
    set(ONOFFhandle,'Position',[40 400 unitwinsize*column unitwinsize*row*2]);
    set(ONOFFhandle,'Name','ONOFF');
else
    figure(ONOFFhandle);
    set(ONOFFhandle,'Position',[40 400 unitwinsize*column unitwinsize*row*2]);
    set(ONOFFhandle,'Name','ONOFF');
end
image(Xpos,Ypos,onoffcolordisp);
axis off;

if newfig
    ONhandle=figure;
    set(ONhandle,'Position',[40 400 unitwinsize*column unitwinsize*row*2]);
    set(ONhandle,'Name','ON');
else
    figure(ONhandle);
    set(ONhandle,'Position',[40 400 unitwinsize*column unitwinsize*row*2]);
    set(ONhandle,'Name','ON');
end
image(Xpos,Ypos,oncolordisp);
axis off;

if newfig
    OFFhandle=figure;
    set(OFFhandle,'Position',[40 400 unitwinsize*column unitwinsize*row*2]);
    set(OFFhandle,'Name','OFF');
else
    figure(OFFhandle);
    set(OFFhandle,'Position',[40 400 unitwinsize*column unitwinsize*row*2]);
    set(OFFhandle,'Name','OFF');
end
image(Xpos,Ypos,offcolordisp);
axis off;

ButtonName=questdlg('Save figures?', 'save option', 'Yes', 'No', 'Yes');
switch ButtonName
    case 'Yes'
    dotpos=strfind(rfpath, '.');
    onsavepath=[rfpath(1:dotpos-1) 'ON' ];
    saveas(ONhandle, onsavepath, 'fig');
    offsavepath=[rfpath(1:dotpos-1) 'OFF' ];
    saveas(OFFhandle, offsavepath, 'fig');
end
