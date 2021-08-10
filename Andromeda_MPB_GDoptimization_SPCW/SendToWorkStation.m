function SendToWorkStation(Workname,Andro_Dir,WorkStation_Dir,cont)
% cont=true;
if cont
%     Andro_Dir=['~/DATA/MPB/20171005_SPCW_param1/'];
%     WorkStation_Dir=['uqml@129.97.41.179:~/Desktop/Youn-MPB-Calculation/20171006_SPCW_candidate2/']
%     Workname='SPCW_ver2'
     
    unix(['scp ' Andro_Dir 'SPCW_ver2.ctl ' WorkStation_Dir 'SPCW_ver2.ctl'])
    unix(['scp ' Andro_Dir 'data.mat ' WorkStation_Dir 'data.mat'])
    unix(['scp ' Andro_Dir 'SPCW_ver2.out ' WorkStation_Dir 'SPCW_ver2.out'])
    unix(['scp ' Andro_Dir '*.h5 ' WorkStation_Dir])
end
