clc;
clear all;

cal_mag_trans_new_ext('result_deltaT.grd',45,10,4,'.\','grd',1);
cal_mag_trans_RTP_ext('result_deltaT.grd',45,10,-45,80,'.\','grd',1);
cal_grad_new_ext('result_deltaT.grd','.\','grd',[0 0 1 0],1);
cal_mag_trans_NSS_ext('result_deltaT.grd',45,10,'.\','grd',1);
cal_mag_trans_RELQ_ext('result_deltaT.grd',45,10,'.\','grd',1);
