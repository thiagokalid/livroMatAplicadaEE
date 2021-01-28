function [id_s0 , iq_s0 , id_r0 , iq_r0] = calc_idq0(VLLrms , Xlr , Xls, Xm , Rr , Rs , s)
% Solução do circuito em regime estacionário, ver exemp 1.1, pag 30
Va = VLLrms * sqrt(2)/ sqrt(3);                         
Zrotor = j*Xlr + Rr/s;                                  
Zm = j*Xm;                                              
Zeq = (Rs + j*Xls) +  (Zm * Zrotor) / (Zm + Zrotor);     
Ia = Va / Zeq;                                        
% Força contra eletromotriz
Ema = Va - (Rs + j*Xls) * Ia;   
% Corrente rotor----Ir
Iraprime = Ema / Zrotor;                               
%*****
Vsp=VLLrms/sqrt(3);
Is=Vsp / Zeq;
Em = Vsp - (Rs + j*Xls) * Is; 
Ir=Em / Zrotor;

% Vectores de tensão e corrente em t=0 no marco de referência do estator
Vs_0 = (3/2) * Va;                            
Is_0 = (3/2) * Ia;                                    
Theta_Is_0 = angle(Is_0);                            
Ir_0 = (-1) * (3/2) * Iraprime;                       
Theta_Ir_0 = angle(Ir_0);

% Alinhamento do eixo a da fase a com o eixo direto d. Assim, Theta_da_0=0
Theta_da_0 = 0;
id_s0 = sqrt(2/3) * abs(Is_0) * cos(Theta_Is_0 - Theta_da_0); 
iq_s0 = sqrt(2/3) * abs(Is_0) * sin(Theta_Is_0 - Theta_da_0); 
id_r0 = sqrt(2/3) * abs(Ir_0) * cos(Theta_Ir_0 - Theta_da_0);
iq_r0 = sqrt(2/3) * abs(Ir_0) * sin(Theta_Ir_0 - Theta_da_0);

end

