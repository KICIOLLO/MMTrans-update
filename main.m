clc;
clear all;

cal_mag_trans_new_ext('result_deltaT.grd',45,10,4,'.\','grd',1);
cal_mag_trans_RTP_ext('result_deltaT.grd',45,10,-45,80,'.\','grd',1);
cal_grad_new_ext('result_deltaT.grd','.\','grd',[0 0 1 0],1);
cal_mag_trans_NSS_ext('result_deltaT.grd',45,10,'.\','grd',1);

cal_magDirection_corr_RTP_NSS('result_deltaT.grd',45,10,1,'.\','grd',1);
cal_magDirection_corr_RTP_NSS('result_deltaT.grd',45,10,2,'.\','grd',1);
cal_magDirection_corr_RTP_NSS('result_deltaT.grd',45,10,3,'.\','grd',1);
cal_magDirection_corr_RTP_NSS('result_deltaT.grd',45,10,4,'.\','grd',1);