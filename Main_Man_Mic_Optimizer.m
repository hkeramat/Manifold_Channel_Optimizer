function Manifold_Opt

Mu=2.8230e-05;
Rho=1.1070;
Cp=1038;
K=0.0422;

Rho_solid=7900;
Cp_solid=500;
K_solid=21;

BC_choice=1;
T_in=300;
FluxInput=10000;
Tbase=400;
material_fin='aluminum';
material_fluid='air';

X_sample=[2000 0.001 2 0.0003 0.003 0.005 300 0.0003 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for q=1
% Calculating fRe and h for single channel model
for i=1:2
    % Reading the variables from X
    % Set the upper and lower limit of Re for linear interpolation of fRe
    Hchn(q,i)=X_sample(q,2);
    alp(q,i)=X_sample(q,3);
    Wchn(q,i)=X_sample(q,4);
    Win(q,i)=X_sample(q,5);
    Hman(q,i)=X_sample(q,6);
    N(q,i)=X_sample(q,7);
    Wman(q,i)=X_sample(q,8);
    Beta(q,i)=X_sample(q,9);

    Win_ratio=Beta(q,i);
    Hin=300e-6;
    Hbase=150e-6;
    Hmnd_ratio=(Hin)/(Hman(q,i)+Hin);
    
    if i==1
            Re_man(q,i)=X_sample(q,1);
        else if last_resid(q,1)>0.01 && Re_man(q,1)>=2000
            Re_man(q,i)=(X_sample(q,1)+Vrb(1,2))/2; 
        else if last_resid(q,1)>0.01 && Re_man(q,1)<2000
            Re_man(q,i)=X_sample(q,1)*2;    
        else if X_sample(q,1)<(Vrb(1,1)+Vrb(1,2))/2
            Re_man(q,i)=X_sample(q,1)+(Vrb(1,1)-Vrb(1,2))/2;
            Rchn=Re_man(q,i)*1/(2*N(q,i))*(Hman(q,i)+Win(q,i))/(Hchn(q,i)+Wchn(q,i));
            Rin=Re_man(q,i)*1/(2*N(q,i))*(Hman(q,i)+Win(q,i))/(Wchn(q,i)+Win(q,i)*Beta(q,i));
            if Rchn<2300 && Rin<2300
            else if Vrb(1,2)==X_sample(q,1)
                    Re_man(q,i)=Vrb(1,2)*2;
                else
                Re_man(q,i)=Vrb(1,2);
                end
            end
            else
            Re_man(q,i)=X_sample(q,1)-(Vrb(1,1)-Vrb(1,2))/2;
            Rchn=Re_man(q,i)*1/(2*N(q,i))*(Hman(q,i)+Win(q,i))/(Hchn(q,i)+Wchn(q,i));
            Rin=Re_man(q,i)*1/(2*N(q,i))*(Hman(q,i)+Win(q,i))/(Wchn(q,i)+Win(q,i)*Beta(q,i));
            if Rchn<2300 && Rin<2300
            else if Vrb(1,2)==X_sample(q,1)
                    Re_man(q,i)=Vrb(1,2)*2;
                else
                Re_man(q,i)=Vrb(1,2);
                end
            end          
            end
        end
        end
   end    

    
    % Channel & manifold property calculation
    t_fin=Wchn(q,i)/alp(q,i);
    Dh_chn(q,i)=2*Hchn(q,i)*Wchn(q,i)/(Hchn(q,i)+Wchn(q,i));
    Dh_man(q,i)=2*Hman(q,i)*Win(q,i)/(Hman(q,i)+Win(q,i));
    m_man(q,i)=Re_man(q,i)*(Hman(q,i)*Win(q,i))*Mu/Dh_man(q,i);
    m_chn(q,i)=m_man(q,i)/(2*N(q,i));
    vc(q,i)=m_chn(q,i)/(Rho*Hchn(q,i)*Wchn(q,i));
    Re_chn(q,i)=Dh_chn(q,i)*m_chn(q,i)/((Hchn(q,i)*Wchn(q,i))*Mu);
    Lflow(q,i)=Wman(q,i)+Win(q,i)-Win(q,i)*Beta(q,i)/2;
    Lchn(q,i)=Wman(q,i)+Win(q,i);
    Lman(q,i)=N(q,i)*(t_fin+Wchn(q,i));
    % Lfd(q,i)=0.05*Re_man(q,i)*Dh_man(q,i);
    if N(q,i)>300
        N_s(q,i)=300;
    else
        N_s(q,i)=N(q,i);
    end
    Lfd(q,i)=N_s(q,i)*(t_fin+Wchn(q,i))/2;
   
    % Set the channel and microchannel dimention to µm for gambit input
    w_chn=Wchn(q,i)*10^6;
    h_chn=Hchn(q,i)*10^6;
    h_base=Hbase*10^6;
    w_in=Win(q,i)*10^6;
    l_man=Hman(q,i)*10^6;
    tfin=t_fin*10^6;
    L_fd=Lfd(q,i)*10^6;
    L_chn=Wman(q,i)*10^6;
    l_man_in=Hin*10^6;
    beta=Beta(q,i);
    
    % Set the grid size 
    if Win_ratio<=0.3
    Grid_Win_fd=16; % 
    Grid_Win_sd=14;
    else if Win_ratio<=0.7
        Grid_Win_fd=19;
        Grid_Win_sd=11;         
        else
            Grid_Win_fd=22;
            Grid_Win_sd=7;         
        end
    end
    Grid_Lchn=21;
    Grid_Hchn=21;
    Grid_Wchn=8;
    Grid_Hbase=5;
    Grid_Hman_top=30;
    Grid_Hman_bottom=8;
    Grid_tfin=6;
    Grid_xfd_in=50;
    Grid_xfd_out=50;
    
    if Win(q,i)>0.006
    if Win_ratio<=0.3
    Grid_Win_fd=16; % 
    Grid_Win_sd=20;
    else if Win_ratio<=0.7
        Grid_Win_fd=22;
        Grid_Win_sd=15;         
        else
            Grid_Win_fd=3;
            Grid_Win_sd=9;         
        end
    end
    Grid_Lchn=31;
    end

    if Win(q,i)>0.012
    if Win_ratio<=0.3
    Grid_Win_fd=18; % 
    Grid_Win_sd=25;
    else if Win_ratio<=0.7
        Grid_Win_fd=24;
        Grid_Win_sd=18;         
        else
            Grid_Win_fd=32;
            Grid_Win_sd=10;         
        end
    end
    Grid_Lchn=42;
    end
    
    if Hman(q,i)>.01
    Grid_Hman_top=40;
    end
    if Hman(q,i)>.02
    Grid_Hman_top=40;
    end    
    
    % Calculate the total Number of Grid 
    % Grid=2*(Grid_Hman*Grid_Win*(2*Grid_tfin+2*Grid_Wchn))+(Grid_Lchn+2*Grid_Win)*(2*Grid_tfin+2*Grid_Wchn)*Grid_Hchn+(Grid_Lchn+2*Grid_Win)*(2*Grid_tfin+2*Grid_Wchn)*Grid_Hbase+(Grid_xfd_out+Grid_xfd_in)*Grid_Hman*Grid_Win;

% Writting Gambit journal file
filename = 'C:/MATLAB5/Journal_Uchn_new_Gam.txt';
gid = fopen(filename, 'w');
fprintf(gid, 'save name "C:/MATLAB5/U_chn_new.dbs" \n');

fprintf(gid, 'vertex create coordinates 0 0 0 \n');
fprintf(gid, 'vertex create coordinates 0 %9.9f 0 \n',l_man);
fprintf(gid, 'vertex create coordinates 0 %9.9f 0 \n',l_man+l_man_in);
fprintf(gid, 'vertex create coordinates 0 %9.9f 0 \n',l_man+l_man_in+h_chn);
fprintf(gid, 'vertex create coordinates 0 %9.9f 0 \n',l_man+l_man_in+h_chn+h_base);
fprintf(gid, 'vertex cmove "vertex.1" "vertex.2" "vertex.3" "vertex.4" "vertex.5" multiple 1 offset %9.9f 0 0 \n',tfin/2);
fprintf(gid, 'vertex cmove "vertex.6" "vertex.7" "vertex.8" "vertex.9" "vertex.10" multiple 1 offset %9.9f 0 0 \n',w_chn/2);
fprintf(gid, 'vertex cmove "vertex.11" "vertex.12" "vertex.13" "vertex.14" "vertex.15" multiple 1 offset %9.9f 0 0 \n',w_chn/2);
fprintf(gid, 'vertex cmove "vertex.16" "vertex.17" "vertex.18" "vertex.19" "vertex.20" multiple 1 offset %9.9f 0 0 \n',tfin/2);
fprintf(gid, 'edge create straight "vertex.21" "vertex.16" "vertex.11" "vertex.6" "vertex.1" "vertex.2" "vertex.3" "vertex.8" "vertex.9" "vertex.14" "vertex.19" "vertex.18" "vertex.23" "vertex.22" \n');
fprintf(gid, 'edge create straight "vertex.22" "vertex.21" \n');
fprintf(gid, 'edge create straight "vertex.3" "vertex.4" "vertex.5" "vertex.10" "vertex.15" "vertex.20" "vertex.25" "vertex.24" \n');
fprintf(gid, 'edge create straight "vertex.23" "vertex.24" \n');
fprintf(gid, 'edge create straight "vertex.8" "vertex.13" "vertex.18" \n');
fprintf(gid, 'edge create straight "vertex.2" "vertex.7" "vertex.12" "vertex.17" "vertex.22" \n');
fprintf(gid, 'edge create straight "vertex.4" "vertex.9" \n');
fprintf(gid, 'edge create straight "vertex.19" "vertex.24" \n');
fprintf(gid, 'face create wireframe "edge.1" "edge.2" "edge.3" "edge.4" "edge.5" "edge.25" "edge.26" "edge.27" "edge.28" "edge.14" real \n');
fprintf(gid, 'face create wireframe "edge.13" "edge.28" "edge.27" "edge.26" "edge.25" "edge.6" "edge.7" "edge.23" "edge.24" "edge.12" real \n');
fprintf(gid, 'face create wireframe "edge.7" "edge.15" "edge.8" "edge.29" real \n');
fprintf(gid, 'face create wireframe "edge.8" "edge.9" "edge.10" "edge.11" "edge.24" "edge.23" real \n');
fprintf(gid, 'face create wireframe "edge.22" "edge.12" "edge.11" "edge.30" real \n');
fprintf(gid, 'face create wireframe "edge.29" "edge.9" "edge.10" "edge.30" "edge.21" "edge.20" "edge.19" "edge.18" "edge.17" "edge.16" real \n');
fprintf(gid, 'volume create translate "face.1" "face.2" "face.4" "face.3" "face.5" "face.6" vector 0 0 %9.9f \n',w_in/2*beta);
fprintf(gid, 'face merge "face.54" "face.55" "face.53" "face.52" \n');
fprintf(gid, 'face merge "face.33" "face.34" \n');
fprintf(gid, 'face merge "face.15" "face.16" \n');
fprintf(gid, 'face merge "face.13" "face.14" \n');
fprintf(gid, 'face merge "face.9" "face.10" \n');
fprintf(gid, 'face merge "face.7" "face.8" \n');
fprintf(gid, 'volume create translate "face.17" vector 0 0 %9.9f \n',w_in/2-w_in/2*beta);
fprintf(gid, 'volume merge "volume.1" "volume.7" real \n');
fprintf(gid, 'face merge "face.63" "face.64" \n');
fprintf(gid, 'face merge "face.65" "face.66" \n');
fprintf(gid, 'face merge "face.57" "face.58" "face.7" \n');
fprintf(gid, 'face merge "face.59" "face.60" "face.9" \n');
fprintf(gid, 'edge create straight "vertex.12" "vertex.11" \n');
fprintf(gid, 'edge create straight "vertex.76" "vertex.71" \n');
fprintf(gid, 'face merge "face.12" "face.62" \n');
fprintf(gid, 'face merge "face.11" "face.61" \n');
fprintf(gid, 'face create wireframe "edge.38" "edge.139" "edge.33" "edge.121" "edge.140" "edge.126" real \n');
fprintf(gid, 'volume split "volume.1" faces "face.68" connected \n');
fprintf(gid, 'volume create translate "face.56" "face.35" "face.40" "face.45" vector 0 0 %9.9f \n',2*(w_in/2-w_in/2*beta)+L_chn);
fprintf(gid, 'volume create translate "face.88" "face.81" "face.98" "face.93" vector 0 0 %9.9f \n',w_in/2*beta);
fprintf(gid, 'face merge "face.78" "face.79" "face.80" "face.77" \n');
fprintf(gid, 'face merge "face.72" "face.73" \n');
fprintf(gid, 'face merge "face.103" "face.104" \n');
fprintf(gid, 'face merge "face.112" "face.113" "face.114" "face.115" \n');
fprintf(gid, 'face merge "face.25" "face.26" \n');
fprintf(gid, 'face merge "face.82" "face.83" \n');
fprintf(gid, 'face merge "face.99" "face.100" \n');
fprintf(gid, 'volume create translate "face.99" "face.122" "face.117" vector 0 %9.9f 0 \n',-l_man_in);
fprintf(gid, 'volume merge "volume.17" "volume.19" "volume.18" real \n');
fprintf(gid, 'face merge "face.129" "face.135" "face.140" \n');
fprintf(gid, 'face merge "face.127" "face.134" "face.139" "face.128" \n');
fprintf(gid, 'face merge "face.132" "face.138" "face.133" "face.143" \n');
fprintf(gid, 'volume create translate "face.129" vector 0 %9.9f 0 \n',-l_man);
fprintf(gid, 'face merge "face.143" "face.144" \n');
fprintf(gid, 'face merge "face.145" "face.146" \n');
fprintf(gid, 'face merge "face.150" "face.151" \n');
fprintf(gid, 'face merge "face.152" "face.153" \n');
fprintf(gid, 'volume create translate "face.143" "face.145" vector 0 0 %9.9f \n',-(w_in/2-w_in/2*beta));
fprintf(gid, 'face merge "face.158" "face.159" \n');
fprintf(gid, 'face merge "face.165" "face.166" \n');
fprintf(gid, 'face merge "face.154" "face.155" \n');
fprintf(gid, 'face merge "face.161" "face.162" \n');
fprintf(gid, 'volume merge "volume.18" "volume.19" "volume.20" real \n');
fprintf(gid, 'edge create straight "vertex.129" "vertex.134" \n');
fprintf(gid, 'edge create straight "vertex.115" "vertex.118" \n');
fprintf(gid, 'face create wireframe "edge.285" "edge.290" "edge.306" "edge.279" "edge.307" "edge.293" real \n');
fprintf(gid, 'volume split "volume.18" faces "face.166" connected \n');
fprintf(gid, 'face merge "face.148" "face.156" \n');
fprintf(gid, 'face merge "face.149" "face.164" \n');
fprintf(gid, 'face merge "face.147" "face.161" \n');
fprintf(gid, 'face merge "face.154" "face.167" \n');
fprintf(gid, 'volume create translate "face.11" vector %9.9f 0 0 \n',-L_fd);
fprintf(gid, 'volume create translate "face.149" vector %9.9f 0 0 \n',L_fd);
fprintf(gid, 'face merge "face.171" "face.174" \n');
fprintf(gid, 'face merge "face.170" "face.173" \n');
fprintf(gid, 'face merge "face.178" "face.181" \n');
fprintf(gid, 'face merge "face.177" "face.180" \n');
fprintf(gid, 'volume merge "volume.3" "volume.10" "volume.13" real \n');
fprintf(gid, 'volume merge "volume.6" "volume.9" "volume.14" real \n');
fprintf(gid, 'volume merge "volume.4" "volume.11" "volume.16" real \n');
fprintf(gid, 'volume merge "volume.5" "volume.12" "volume.15" real \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge picklink "edge.199" "edge.198" "edge.176" "edge.175" "edge.80" "edge.79"  "edge.234" "edge.228" "edge.189" "edge.182" "edge.97" "edge.88" "edge.22" "edge.15" "edge.11" "edge.8" \n');
fprintf(gid, 'edge mesh "edge.8" "edge.11" "edge.15" "edge.22" "edge.79" "edge.80" "edge.88" "edge.97" "edge.175" "edge.176" "edge.182" "edge.189" "edge.198" "edge.199" "edge.228" "edge.234" successive ratio1 1.2 ratio2 1.2 intervals %9.9f \n',Grid_Hchn);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge picklink "edge.217" "edge.216" "edge.161" "edge.162" "edge.114" "edge.113" "edge.21" "edge.16" \n');
fprintf(gid, 'edge mesh "edge.16" "edge.21" "edge.113" "edge.114" "edge.162" "edge.161" "edge.216" "edge.217" successive ratio1 1 intervals %9.9f \n',Grid_Hbase);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge picklink "edge.261" "edge.253" "edge.256" "edge.248" "edge.66" "edge.65" "edge.13" "edge.6" \n');
fprintf(gid, 'edge mesh "edge.6" "edge.13" "edge.65" "edge.66" "edge.248" "edge.256" "edge.253" "edge.261" successive ratio1 1.2 ratio2 1.2 intervals %9.9f \n',Grid_Hman_bottom);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge picklink "edge.285" "edge.336" "edge.327" "edge.281" "edge.279" "edge.277" "edge.297" "edge.134" "edge.284" "edge.146" "edge.133" "edge.14" "edge.145" "edge.5" "edge.324" "edge.315" \n');
fprintf(gid, 'edge mesh "edge.315" "edge.324" "edge.5" "edge.145" "edge.14" "edge.133" "edge.146" "edge.284" "edge.134" "edge.285" "edge.297" "edge.277" "edge.279" "edge.281" "edge.327" "edge.336" successive ratio1 1.2 ratio2 1.2 intervals %9.9f \n',Grid_Hman_top);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge picklink "edge.334" "edge.333" "edge.313" "edge.312" "edge.272" "edge.271" "edge.258" "edge.250" "edge.229" "edge.223" "edge.211" "edge.207" "edge.206" "edge.202" "edge.195" "edge.193" "edge.192" "edge.190" "edge.318" "edge.317" "edge.108" "edge.104" "edge.94" "edge.85" "edge.76" "edge.74" "edge.60" "edge.59" "edge.57" "edge.56" "edge.40" "edge.38" "edge.36" "edge.35" "edge.33" "edge.31" \n');
fprintf(gid, 'edge mesh "edge.31" "edge.33" "edge.35" "edge.36" "edge.38" "edge.40" "edge.56" "edge.57" "edge.59" "edge.60" "edge.74" "edge.76" "edge.85" "edge.94" "edge.104" "edge.108" "edge.317" "edge.318" "edge.190" "edge.192" "edge.193" "edge.195" "edge.202" "edge.206" "edge.207" "edge.211" "edge.223" "edge.229" "edge.250" "edge.258" "edge.271" "edge.272" "edge.312" "edge.313" "edge.333" "edge.334" successive ratio1 1 intervals %9.9f \n',Grid_Win_fd);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge picklink "edge.330" "edge.329" "edge.305" "edge.302" "edge.293" "edge.291" "edge.290" "edge.288" "edge.322" "edge.321" "edge.128" "edge.126" "edge.124" "edge.123" "edge.121" "edge.119" \n');
fprintf(gid, 'edge mesh "edge.119" "edge.121" "edge.123" "edge.124" "edge.126" "edge.128" "edge.321" "edge.322" "edge.288" "edge.290" "edge.291" "edge.293" "edge.302" "edge.305" "edge.329" "edge.330" successive ratio1 1.2 ratio2 1.2 intervals %9.9f \n',Grid_Win_sd);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge picklink "edge.184" "edge.177" "edge.169" "edge.167" "edge.156" "edge.152" "edge.151" "edge.150" "edge.148" "edge.147" \n');
fprintf(gid, 'edge mesh "edge.147" "edge.148" "edge.150" "edge.151" "edge.152" "edge.156" "edge.167" "edge.169" "edge.177" "edge.184" successive ratio1 1.2 ratio2 1.2 intervals %9.9f \n',Grid_Lchn);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge picklink "edge.49" "edge.48" "edge.201" "edge.200" "edge.159" "edge.158" "edge.82" "edge.81" "edge.69" "edge.68" "edge.294" "edge.275" "edge.131" "edge.2" "edge.283" "edge.274" "edge.130" "edge.3" "edge.298" "edge.243" "edge.236" "edge.137" "edge.27" "edge.287" "edge.242" "edge.235" "edge.136" "edge.26" "edge.220" "edge.219" "edge.197" "edge.196" "edge.174" "edge.173" "edge.165" "edge.164" "edge.117" "edge.116" "edge.24" "edge.23" "edge.19" "edge.18" "edge.10" "edge.9" \n');
fprintf(gid, 'edge mesh "edge.9" "edge.10" "edge.18" "edge.19" "edge.23" "edge.24" "edge.68" "edge.69" "edge.81" "edge.82" "edge.116" "edge.117" "edge.158" "edge.159" "edge.164" "edge.165" "edge.173" "edge.174" "edge.196" "edge.197" "edge.200" "edge.201" "edge.219" "edge.220" "edge.26" "edge.48" "edge.136" "edge.235" "edge.242" "edge.287" "edge.27" "edge.49" "edge.137" "edge.236" "edge.243" "edge.298" "edge.3" "edge.130" "edge.274" "edge.283" "edge.2" "edge.131" "edge.275" "edge.294" successive ratio1 1 intervals %9.9f \n',Grid_Wchn);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge delete "edge.4" "edge.129" "edge.273" "edge.282" "edge.1" "edge.132" "edge.276" "edge.295" "edge.25" "edge.47" "edge.135" "edge.247" "edge.252" "edge.286" "edge.28" "edge.50" "edge.138" "edge.254" "edge.259" "edge.299" "edge.7" "edge.17" "edge.29" "edge.67" "edge.90" "edge.115" "edge.157" "edge.163" "edge.181" "edge.212" "edge.218" "edge.233" "edge.12" "edge.20" "edge.30" "edge.70" "edge.98" "edge.118" "edge.160" "edge.166" "edge.187" "edge.215" "edge.221" "edge.226" keepsettings onlymesh \n');
fprintf(gid, 'edge picklink "edge.215" "edge.160" "edge.98" "edge.70" "edge.212" "edge.157" "edge.90" "edge.67" "edge.50" "edge.47" "edge.226" "edge.221" "edge.187" "edge.166" "edge.118" "edge.30" "edge.20" "edge.12" "edge.233" "edge.218" "edge.181" "edge.163" "edge.115" "edge.29" "edge.17" "edge.7" "edge.299" "edge.259" "edge.254" "edge.138" "edge.28" "edge.286" "edge.252" "edge.247" "edge.135" "edge.25" "edge.295" "edge.276" "edge.132" "edge.1" "edge.282" "edge.273" "edge.129" "edge.4" \n');
fprintf(gid, 'edge mesh "edge.4" "edge.129" "edge.273" "edge.282" "edge.1" "edge.132" "edge.276" "edge.295" "edge.25" "edge.47" "edge.135" "edge.247" "edge.252" "edge.286" "edge.28" "edge.50" "edge.138" "edge.254" "edge.259" "edge.299" "edge.7" "edge.17" "edge.29" "edge.67" "edge.90" "edge.115" "edge.157" "edge.163" "edge.181" "edge.212" "edge.218" "edge.233" "edge.12" "edge.20" "edge.30" "edge.70" "edge.98" "edge.118" "edge.160" "edge.166" "edge.187" "edge.215" "edge.221" "edge.226" successive ratio1 1 intervals %9.9f \n',Grid_tfin);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge delete "edge.314" "edge.316" "edge.323" "edge.325" keepsettings onlymesh \n');
fprintf(gid, 'edge picklink "edge.325" "edge.323" "edge.316" "edge.314" \n');
fprintf(gid, 'edge mesh "edge.314" "edge.316" "edge.323" "edge.325" successive ratio1 1.2 ratio2 1.2 intervals %9.9f \n',Grid_xfd_in);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'undo begingroup \n');
fprintf(gid, 'edge delete "edge.326" "edge.328" "edge.335" "edge.337" keepsettings onlymesh \n');
fprintf(gid, 'edge picklink "edge.337" "edge.335" "edge.328" "edge.326" \n');
fprintf(gid, 'edge mesh "edge.326" "edge.328" "edge.335" "edge.337" successive ratio1 1.2 ratio2 1.2 intervals %9.9f \n',Grid_xfd_out);
fprintf(gid, 'undo endgroup \n');
fprintf(gid, 'face mesh "face.1" "face.2" "face.3" "face.4" "face.5" "face.6" "face.11" "face.12" "face.13" "face.15" "face.22" "face.23" "face.24" "face.25" "face.27" "face.28" "face.31" "face.32" "face.33" "face.37" "face.39" "face.43" "face.44" "face.50" "face.51" "face.54" "face.57" "face.59" "face.63" "face.65" "face.67" "face.68" "face.69" "face.70" "face.71" "face.72" "face.74" "face.75" "face.76" "face.78" "face.82" "face.84" "face.85" "face.89" "face.90" "face.94" "face.96" "face.99" "face.101" "face.102" "face.103" "face.105" "face.106" "face.109" "face.110" "face.111" "face.112" "face.116" "face.117" "face.119" "face.121" "face.122" "face.123" "face.126" "face.127" "face.129" "face.132" "face.136" "face.142" "face.147" "face.148" "face.149" "face.150" "face.152" "face.153" "face.154" "face.158" "face.160" "face.165" "face.166" "face.168" "face.169" "face.170" "face.171" "face.172" "face.175" "face.176" "face.177" "face.178" "face.179" "face.182" map size 1 \n');
fprintf(gid, 'volume mesh "volume.2" "volume.1" "volume.8" "volume.17" "volume.18" "volume.19" "volume.20" "volume.21" "volume.3" "volume.6" "volume.4" "volume.5" map size 1 \n');
fprintf(gid, 'volume merge "volume.6" "volume.4" "volume.5" real \n');
fprintf(gid, 'volume mesh "volume.6" submap size 1 \n');
fprintf(gid, 'physics create "Microchannel_in" btype "INTERIOR" face "face.99" \n');
fprintf(gid, 'physics create "Microchannel_out" btype "INTERIOR" face "face.25" \n');
fprintf(gid, 'physics create "CV_center_in" btype "INTERIOR" face "face.166" \n');
fprintf(gid, 'physics create "CV_center_out" btype "INTERIOR" face "face.68" \n');
fprintf(gid, 'physics create "Inlet" btype "MASS_FLOW_INLET" face "face.179" \n');
fprintf(gid, 'physics create "Outlet" btype "PRESSURE_OUTLET" face "face.172" \n');
fprintf(gid, 'physics create "Flux" btype "WALL" face "face.78" "face.112" "face.54" \n');
fprintf(gid, 'physics create "Wall" btype "WALL" face "face.85" "face.102" "face.32" "face.84" "face.101" "face.31" "face.72" "face.103" "face.33" \n');
fprintf(gid, 'physics create "Symmetry" btype "SYMMETRY" face "face.150" "face.152" "face.182" "face.132" "face.116" "face.126" "face.105" "face.121" "face.4" "face.6" "face.3" "face.5" "face.2" "face.1" "face.69" "face.169" \n');
fprintf(gid, 'physics modify "Symmetry" btype face "face.150" "face.152" "face.182" "face.132" "face.116" "face.126" "face.105" "face.121" "face.4" "face.6" "face.3" "face.5" "face.2" "face.1" "face.69" "face.169" \n');
fprintf(gid, 'physics create "Interior" btype "INTERIOR" face "face.149" "face.11" "face.129" "face.168" "face.13" "face.68" "face.15" \n');
fprintf(gid, 'physics create "water" ctype "FLUID" volume "volume.2" "volume.1" "volume.8" "volume.17" "volume.18" "volume.19" "volume.20" "volume.21" "volume.3" \n');
fprintf(gid, 'physics create "solid" ctype "SOLID" volume "volume.6" \n');

fprintf(gid, 'export fluent5 "C:/MATLAB5/U_chn_new.msh" \n');
fprintf(gid, 'save \n');
fclose(gid);

pause(2)

% Running Gambit (Creating mesh file)
! gambit -inp "C:/MATLAB5/Journal_Uchn_new_Gam.txt"

pause(2)

%Tbase
% Printing Fluent journal file
Scale=0.000001;
Iter=300;
mchn=m_chn(q,i);
filename = 'C:/MATLAB5/Journal_Uchn_new_Fln.txt';
fid = fopen(filename, 'w');
fprintf(fid, 'file read-case C:/MATLAB5/U_chn_new.msh \n');
fprintf(fid, 'grid/scale/ %1.6f %1.6f %1.6f \n',Scale,Scale,Scale);
fprintf(fid, 'define/models/energy yes no no no yes \n');
fluid_check=strcmp(material_fluid,'air');
if fluid_check==1
fprintf(fid, 'define/boundary-conditions/fluid water no no no yes -1 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 1 none no no no \n');
fprintf(fid, 'q \n');
else
fprintf(fid, 'define/materials/copy fluid %0.18s \n',material_fluid);
fprintf(fid, 'define/boundary-conditions/fluid water yes %0.18s no no yes -1 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 1 none no no no \n',material_fluid);
fprintf(fid, 'q \n');
end
solid_check=strcmp(material_fin,'inconel');
if solid_check==1
fprintf(fid, 'define/materials/ change-create aluminum inconel yes constant %6.9f yes constant %6.9f yes constant %6.9f yes \n',Rho_solid,Cp_solid,K_solid);    
fprintf(fid, 'define/boundary-conditions/solid solid no no no yes -1 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 1 none no no no \n');
fprintf(fid, 'q \n');
else
fprintf(fid, 'define/materials/copy solid %0.9s \n',material_fin);
fprintf(fid, 'define/boundary-conditions/solid solid yes %0.9s no no yes -1 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 1 none no no no \n',material_fin);
fprintf(fid, 'q \n');
end
fprintf(fid, 'define/material/change-create/air/air yes constant %6.9f yes constant %6.9f yes constant %6.9f yes constant %6.9f no no no \n',Rho,Cp,K,Mu);

fprintf(fid, 'define/boundary-conditions/mass-flow inlet yes yes no %1.10f no %4.1f no 0 no yes \n',mchn,T_in); 
fprintf(fid, 'q \n');
fprintf(fid, 'define/boundary-conditions/pressure-outlet outlet no 0 no 300 no yes no no no \n');
fprintf(fid, 'q \n');
if BC_choice==1
fprintf(fid, 'define/boundary-conditions/wall flux 0 no 0 yes %0.9s no no %8.1f no no 1 \n',material_fin,FluxInput); 
fprintf(fid, 'q \n');
else if BC_choice==0
fprintf(fid, 'define/boundary-conditions/wall flux 0 no 0 yes %0.9s yes temperature no %8.1f no no 1 \n',material_fin,Tbase);
fprintf(fid, 'q \n');
end
end
fprintf(fid, 'solve/monitors/residual/convergence-criteria 0.00001 0.00001 0.00001 0.00001 0.000001  \n');
fprintf(fid, 'solve/initialize/hyb-initialization \n');

fprintf(fid, 'define/models/energy no \n');
fprintf(fid, 'solve/iterate %3.0f \n',3);
fprintf(fid, 'solve/iterate %3.0f \n',350);
fprintf(fid, 'define/models/energy yes no no no yes \n');

fprintf(fid, 'solve/monitors/surface/set-monitor Pressure_result1 "Area-Weighted Average" pressure Microchannel_in () no yes yes C:/MATLAB5/Pressure_result1b 20 \n');
fprintf(fid, 'solve/monitors/surface/set-monitor Pressure_result2 "Area-Weighted Average" pressure Microchannel_out () no yes yes C:/MATLAB5/Pressure_result2b 20 \n');
fprintf(fid, 'solve/monitors/surface/set-monitor Pressure_result3 "Area-Weighted Average" pressure CV_center_in () no yes yes C:/MATLAB5/Pressure_result3b 20 \n');
fprintf(fid, 'solve/monitors/surface/set-monitor Pressure_result4 "Area-Weighted Average" pressure CV_center_out () no yes yes C:/MATLAB5/Pressure_result4b 20 \n');
fprintf(fid, 'solve/monitors/surface/set-monitor Inlet_T "Area-Weighted Average" temperature inlet () no yes yes C:/MATLAB5/Inlet_Tb 20 \n');
fprintf(fid, 'solve/monitors/surface/set-monitor Base_T "Area-Weighted Average" temperature flux () no yes yes C:/MATLAB5/Base_Tb 20 \n');
fprintf(fid, 'solve/monitors/surface/set-monitor Wall_T "Area-Weighted Average" temperature wall () no yes yes C:/MATLAB5/Wall_Tb 20 \n');
fprintf(fid, 'solve/monitors/surface/set-monitor flux "Area-Weighted Average" heat-flux flux () no yes yes C:/MATLAB5/fluxb 20 \n');
fprintf(fid, 'solve/monitors/surface/set-monitor flux2 "Area-Weighted Average" heat-flux wall () no yes yes C:/MATLAB5/flux2b 20 \n');

%fprintf(fid, 'solve/monitors/residual/normalize? yes \n');
fprintf(fid, 'solve/monitors/residual/plot? yes \n');
fprintf(fid, 'solve/iterate 10 \n');
fprintf(fid, 'solve/iterate 150 \n');
fprintf(fid, 'file write-case-data C:/MATLAB5/U_chn_new1 yes \n');
fprintf(fid, 'plot/residuals-set plot-to-file C:/MATLAB5/Resid1 \n');
fprintf(fid, 'plot/residuals yes yes yes yes yes \n');
fprintf(fid, 'plot/residuals-set plot-to-file C:/MATLAB5/Signal \n \n');
fprintf(fid, 'exit \n');
fclose(fid);

pause(2)


% Running Fluent (Creating .cas and .data file)
! fluent 3ddp -i "C:/MATLAB5/Journal_Uchn_new_Fln.txt"


% Ask mathlab to wait until fluent running finish
for iterate=1:1000000
try              
[Signal] = textread('C:/MATLAB5/Signal','%s', 1000);  
catch error  
  pause(1)
continue       
end
pause(2)
break
end

pause(8)

% Check Confergance
[conf1,conf2] = textread('C:\MATLAB5\Resid1','%s %s', 1000);
confergance1=str2double(conf1);
confergance2=str2double(conf2);
conf_check=max(confergance1);
last_resid(q,i)=confergance2(conf_check)

% Check Tbase 
[iter1,Temperature2] = textread('C:\MATLAB5\Base_Tb','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Temperature_2=str2double(Temperature2);
T_base_avg=Temperature_2(iter_max)

if last_resid(q,1)>0.01 || T_base_avg>1000
    last_resid(q,2)=10;
    delete('C:/MATLAB5/Journal_Uchn_new_Gam.txt')
    delete('C:/MATLAB5/U_chn_new.msh')
    delete('C:/MATLAB5/Journal_Uchn_new_Fln.txt')
    delete('C:/MATLAB5/Pressure_result1b')
    delete('C:/MATLAB5/Pressure_result2b')
    delete('C:/MATLAB5/Pressure_result3b')
    delete('C:/MATLAB5/Pressure_result4b')
    delete('C:/MATLAB5/Pressure_result5b')
    delete('C:/MATLAB5/Pressure_result6b')
    delete('C:/MATLAB5/Pressure_result7b')
    delete('C:/MATLAB5/Pressure_result8b')
    delete('C:/MATLAB5/Inlet_Tb')
    delete('C:/MATLAB5/Wall_Tb')
    delete('C:/MATLAB5/base_Tb')
    delete('C:/MATLAB5/fluxb')
    delete('C:/MATLAB5/flux2b')
    delete('C:/MATLAB5/Signal')
    delete('C:/MATLAB5/Resid1')
    delete('C:/MATLAB5/default_id.lok')
    delete('C:/Matlab Interface/default_id.lok')
    break
else

% Post processing reading the Fluent simulation results (For uniform mesh)
% Calculating pressure at microchannel channel inlet
[iter1,Pressure1] = textread('C:\MATLAB5\Pressure_result1b','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Pressure_in=str2double(Pressure1);
P_mchn_in_avg(q,i)=Pressure_in(iter_max)

% Calculating pressure at microchannel channel outlet
[iter1,Pressure2] = textread('C:\MATLAB5\Pressure_result2b','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Pressure_out=str2double(Pressure2);
P_mchn_out_avg(q,i)=Pressure_out(iter_max);

% Calculating pressure at CV inlet middle
[iter1,Pressure3] = textread('C:\MATLAB5\Pressure_result3b','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Pressure_in2=str2double(Pressure3);
P_CV_in_avg(q,i)=Pressure_in2(iter_max)

% Calculating pressure at CV outlet middle
[iter1,Pressure4] = textread('C:\MATLAB5\Pressure_result4b','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Pressure_out2=str2double(Pressure4);
P_CV_out_avg(q,i)=Pressure_out2(iter_max);

% Calculating pressure drop at microchannel and CV_middle 
dPmchn(q,i)=P_mchn_in_avg(q,i)-P_mchn_out_avg(q,i);
dP12(q,i)=P_CV_in_avg(q,i)-P_CV_out_avg(q,i)

% Calculating fRe at microchannel and CV_middle
fRe12(q,i)=Re_chn(q,i)*dP12(q,i)*Dh_chn(q,i)/(4*Lflow(q,i))*2/(Rho*vc(q,i)^2);
fRe_chn(q,i)=Re_chn(q,i)*dPmchn(q,i)*Dh_chn(q,i)/(4*Lflow(q,i))*2/(Rho*vc(q,i)^2);

% Calculating P/V (pumping power/Volume) at microchannel
Volume(q,i)=(t_fin+Wchn(q,i))*Lchn(q,i)*(Hchn(q,i)+Hbase+Hman(q,i)+Hin);
P_pump(q,i)=m_chn(q,i)*dPmchn(q,i)/Rho/Volume(q,i);

% Calculating Tin 
[iter1,Temperature1] = textread('C:\MATLAB5\Inlet_Tb','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Temperature_1=str2double(Temperature1);
T_in_avg=Temperature_1(iter_max);

% Calculating Tbase 
[iter1,Temperature2] = textread('C:\MATLAB5\Base_Tb','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Temperature_2=str2double(Temperature2);
T_base_avg=Temperature_2(iter_max)

% Calculating Twall 
[iter1,Temperature3] = textread('C:\MATLAB5\Wall_Tb','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Temperature_3=str2double(Temperature3);
T_wall_avg=Temperature_3(iter_max);

% Calculating Heat flux at the base 
[iter1,Heat] = textread('C:\MATLAB5\fluxb','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Heat_tot=str2double(Heat);
Q_tot_base(q,i)=Heat_tot(iter_max);

% Calculating Heat flux at the wall 
[iter1,Heat2] = textread('C:\MATLAB5\flux2b','%s %s', 1000);
iter=str2double(iter1);
iter_max0=size(iter);
iter_max=iter_max0(1,1);
Heat_tot_2=str2double(Heat2);
Q_tot_wall(q,i)=Heat_tot_2(iter_max);

% Calculating heat transfer coefficient at wall and base (h) and heat flow
% rate (Q)
h_wall(q,i)=Q_tot_wall(q,i)/(T_wall_avg-T_in_avg);
h_Base1(q,i)=Q_tot_base(q,i)/(T_base_avg-T_in_avg);
h_BaseV(q,i)=Q_tot_base(q,i)/(T_base_avg-T_in_avg)/(Hchn(q,i)+Hbase+Hman(q,i)+Hin)
Heat_flux(q,i)=Q_tot_base(q,i)*(t_fin+Wchn(q,i))*Lchn(q,i);

%pause(1)

% Deleting gambit and fluent files
delete('C:/MATLAB5/Journal_Uchn_new_Gam.txt')
delete('C:/MATLAB5/U_chn_new.msh')
delete('C:/MATLAB5/Journal_Uchn_new_Fln.txt')
delete('C:/MATLAB5/Pressure_result1b')
delete('C:/MATLAB5/Pressure_result2b')
delete('C:/MATLAB5/Pressure_result3b')
delete('C:/MATLAB5/Pressure_result4b')
delete('C:/MATLAB5/Pressure_result5b')
delete('C:/MATLAB5/Pressure_result6b')
delete('C:/MATLAB5/Pressure_result7b')
delete('C:/MATLAB5/Pressure_result8b')
delete('C:/MATLAB5/Inlet_Tb')
delete('C:/MATLAB5/Wall_Tb')
delete('C:/MATLAB5/base_Tb')
delete('C:/MATLAB5/fluxb')
delete('C:/MATLAB5/flux2b')
delete('C:/MATLAB5/Signal')
delete('C:/MATLAB5/Resid1')
delete('C:/MATLAB5/default_id.lok')
delete('C:/Matlab Interface/default_id.lok')
save('C:/Matlab Interface/Trial_1')
end
end
end
%! fluent 3ddp

for q=1
q
% Calculating fRe and h for single channel model
if last_resid(q,1)>0.01 || last_resid(q,2)>0.01 || last_resid(q,1)==0 || last_resid(q,2)==0
    for iteration=1:2
    z1=z1+1;
    P_mchn_in_avg1(z1)=0;
    P_mchn_out_avg1(z1)=0;
    P_CV_in_avg1(z1)=0;
    P_CV_out_avg1(z1)=0;
    Q_tot_base1(z1)=0;
    Q_tot_wall1(z1)=0;
    fRe12V1(z1)=0;
    fRe_chnV1(z1)=0;
    vc1(z1)=0;
    a_fRe1(z1)=0;
    b_fRe1(z1)=0;
    a_hb1(z1)=0;
    b_hb1(z1)=0;
    a_f12_1(z1)=0;
    b_f12_1(z1)=0;

    h_Base_1(z1)=0;
    h_wall1(z1)=0;
    P_pump1(z1)=0;
    Heat_flux1(z1)=0;
    Volume_tot1(z1)=0;
    dPt_num1(z1)=0;
    F1_mc1(z1)=0;
    F2_mc1(z1)=0;
    F1_vc1(z1)=0;
    F2_vc1(z1)=0;
    F1_hc1(z1)=0;
    F2_hc1(z1)=0;
    F1_dPmch1(z1)=0;
    F2_dPmch1(z1)=0;

    Re_man1(z1)=0;
    Hchn1(z1)=0;
    alp1(z1)=0;
    Wchn1(z1)=0;
    Win1(z1)=0;
    Hman1(z1)=0;
    N1(z1)=0;
    Wman1(z1)=0;
    P_pump_tot1(z1)=0;
    h_BaseV1(z1)=0;
    COP1(z1)=0;
     if SelectFinType==2
         beta1(z1)=0;
         gamma1(z1)=0;
     end
    end
    Bad_data_count1=Bad_data_count1+1;
else

% Start analytical solver code for full model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating linear eq for fRe12 ==> fRe12=a*vc+b
b=(fRe12(q,2)-fRe12(q,1)*(vc(q,2)/vc(q,1)))/(1-(vc(q,2)/vc(q,1)));
a=(fRe12(q,1)-b)/vc(q,1);

% Calculating linear eq for h_base ==> h_base=a1*vc+b1
b1=(h_Base1(q,2)-h_Base1(q,1)*(vc(q,2)/vc(q,1)))/(1-(vc(q,2)/vc(q,1)));
a1=(h_Base1(q,1)-b1)/vc(q,1);

% Calculating linear eq for fRe_chn ==> fRe_chn=a2*vc+b2
b2=(fRe_chn(q,2)-fRe_chn(q,1)*(vc(q,2)/vc(q,1)))/(1-(vc(q,2)/vc(q,1)));
a2=(fRe_chn(q,1)-b2)/vc(q,1);

% Loop to calculate dP_tot for lower and upper boundary of Re
for i=1:2
% Update z    
z1=z1+1;

% The differntial equation is on the form of:
% d2v1/dx2(k3*dv1/dx-1)+2k1*dv1/dx+2k1k2v1-k1k2=0
% BC: v1(x=0)=1 & v1(x=1)=0
% Where:
% k1=((2*N(q,i)*Ac*Rho*v_in)/Ah)*(Dh_chn(q,i)^2/(2*Mu*Lflow(q,i)*b));
% k2=Ph*fRe_h*Lman(q,i)*Mu/(2*Ah*Rho*Dh_man(q,i)*v_in);
% k3=(a/b)*(2*Ah*v_in)/(2*N(q,i)*Ac);
% Note: v1 and x is normilize velocity and distance wrt vin and Lchn

Ah=Win(q,i)*Hman(q,i);
Ac=Hchn(q,i)*Wchn(q,i);
Ph=2*(Win(q,i)+Hman(q,i));
v_in=m_man(q,i)/(Rho*Ah);

% Calculating horizontal pressure drop at manifold using fully develop
% Corrolation
a_max=max(Win(q,i),Hman(q,i));
a_min=min(Win(q,i),Hman(q,i));
alp_h=a_min/a_max;
fRe_h=24*(1-1.3553*alp_h+1.9467*alp_h^2-1.7012*alp_h^3+0.9564*alp_h^4-0.2537*alp_h^5);

% Calculating k1, k2, and k3
k1=((2*N(q,i)*Ac*Rho*v_in)/Ah)*(Dh_chn(q,i)^2/(2*Mu*Lflow(q,i)*b));
k2=Ph*Lman(q,i)*Mu/(2*Ah*Rho*Dh_man(q,i)*v_in);
k3=(a/b)*(2*Ah*v_in)/(2*N(q,i)*Ac);
Beta1=0.8;
Beta2=-0.8;
fRe_h=24*(1-1.3553*alp_h+1.9467*alp_h^2-1.7012*alp_h^3+0.9564*alp_h^4-0.2537*alp_h^5);

% Solving the differential equation numerically using bvp4c for v1
nx=1000; nxp=nx+1; L=1; dx=L/nx;
xrange=linspace(0,L,nxp);
guess=[.5, .5];
struct1=bvpinit(xrange, guess);
bvsol=bvp4c(@f_tur, @BC_tur, struct1, [],k1,k2,k3,Rho,Mu,Dh_man(q,i),v_in,Beta1,Beta2,fRe_h);
v1_num=deval(bvsol, xrange);
plot(xrange, v1_num(1,:))

% Solving dp/dx using dp/dx=-2*v1*dv1/dx-k2*v1
for j=1:nx
x2(j)=j/nx;

Re1=Rho*v1_num(1,j)*Dh_man(q,i)*v_in/Mu;
if Re1>2300
fRe1=Re1*(0.79*log(Re1)-1.64)^-2/4;
else
    fRe1=fRe_h;
end

dPdx_1(j)=-2*v1_num(1,j)*(v1_num(1,j+1)-v1_num(1,j))/dx-fRe1*k2*v1_num(1,j)+Beta1*v1_num(1,j)*(v1_num(1,j+1)-v1_num(1,j))/dx;
end

% Solving vc at x=1 ==> vc=-Ah*vin/(2*N*Ac)*dv1/dx
vc3(z1)=-Ah/(2*Ac)*1/N(q,i)*(v1_num(1,nxp)-v1_num(1,nx))/(L/nx)*v_in;
% Solving p_mnd(x=1) using integration of dp/dx 
% Note dp/dx is for normzlize P so it need to be multiple by Rho*v_in^2
P1_1(z1)=trapz(x2,dPdx_1)*Rho*v_in^2;
% dp12(x=1) calculatuion
dP12_1(z1)=1/k1*(1/2*k3*((v1_num(1,nxp)-v1_num(1,nx))/(dx))^2-((v1_num(1,nxp)-v1_num(1,nx))/(dx)))*Rho*v_in^2;
% dpchn(x=1) calculation
dPchn(z1)=2*Mu*Lflow(q,i)/Dh_chn(q,i)^2*(a2*vc3(z1)^2+b2*vc3(z1));
dP12_2(z1)=2*Mu*Lflow(q,i)/Dh_chn(q,i)^2*(a*vc3(z1)^2+b*vc3(z1));
% dptot calculation ==> dptot=dp12(x=1)-dp_mnd(x=1)
% Note: for the analytical solution, it is assumed p=0 at inlet and
% increase to outlet ==> dp_mnd=p_mnd 
% (Note: here + dP meant pressure increase)
dPt_num1(z1)=-(P1_1(z1)-dP12_1(z1));
% Calculate the ratio pressure drop in mirochannel and manifold
Ratio(z1)=dPchn(z1)/(dPt_num1(z1)-dPchn(z1));
% Calculate total pumping power
Volume_tot(z1)=Volume(q,i)*N(q,i);
%mass(i)=N(q,i)*t_fin*Lchn(q,i)*Hchn(q,i)+N(q,i)*(t_fin+Wchn(q,i))*Lchn(q,i)*Hbase+N(q,i)*(t_fin+Wchn(q,i))*Wman(q,i)*(Hman(q,i))+N(q,i)*(t_fin+Wchn(q,i))*Hin*(Lchn(q,i)-Beta(q,i)*Win(q,i));
P_pump_tot1(z1)=m_man(q,i)/2*dPt_num1(z1)/Rho/Volume_tot(z1)

% Saving the result for z
P_mchn_in_avg1(z1)=P_mchn_in_avg(q,i);
P_mchn_out_avg1(z1)=P_mchn_out_avg(q,i);
P_CV_in_avg1(z1)=P_CV_in_avg(q,i);
P_CV_out_avg1(z1)=P_CV_out_avg(q,i);
Q_tot_base1(z1)=Q_tot_base(q,i);
Q_tot_wall1(z1)=Q_tot_wall(q,i);
fRe12V1(z1)=fRe12(q,i);
fRe_chnV1(z1)=fRe_chn(q,i);
vc1(z1)=vc(q,i);
a_fRe1(z1)=a;
b_fRe1(z1)=b;
a_hb1(z1)=a1;
b_hb1(z1)=b1;
a_f12_1(z1)=a2;
b_f12_1(z1)=b2;

Re_man1(z1)=Re_man(q,i);
Hchn1(z1)=Hchn(q,i);
alp1(z1)=alp(q,i);
Wchn1(z1)=Wchn(q,i);
Win1(z1)=Win(q,i);
Hman1(z1)=Hman(q,i);
N1(z1)=N(q,i);
Wman1(z1)=Wman(q,i);
beta1(z1)=Beta(q,i);
h_Base_1(z1)=h_Base1(q,i);
h_BaseV1(z1)=h_BaseV(q,i)
h_wall1(z1)=h_wall(q,i);
P_pump1(z1)=P_pump(q,i);
Heat_flux1(z1)=Heat_flux(q,i);
Volume_tot1(z1)=Volume_tot(z1);
COP1(z1)=h_BaseV1(z1)/P_pump_tot1(z1);


% Calculating the variation of vc, hc, and dPmchn using F1 and F2 where:
% F1=(Max-min)/Max
% F2=Standard_deviation/Mean
for k=1:nx
    vc_analytical(k)=(-Ah/(2*Ac)*1/N(q,i)*(v1_num(1,k+1)-v1_num(1,k))/(L/nx))*v_in;
    hc_analytical(k)=a1*vc_analytical(k)+b1;
    dPmch_analytical(k)=a2*vc_analytical(k)+b2;
    mc(k)=vc_analytical(k)*Rho*Ac;
    xprime(k)=(k)/nx*Lman(q,i);
end

%h_Base_analytical(z1)=mean(hc_analytical)

for k=1:nx
    f2_mc(k)=(mc(k)-mean(mc))^2;
    f2_vc(k)=(vc_analytical(k)-mean(vc_analytical))^2;
    f2_dPmch(k)=(dPmch_analytical(k)-mean(dPmch_analytical))^2;
    f2_hc(k)=(hc_analytical(k)-mean(hc_analytical))^2;
end

F1_mc1(z1)=(max(mc)-min(mc))/max(mc);
F2_mc1(z1)=sqrt(sum(f2_mc)/nx)/(mean(mc));

F1_vc1(z1)=(max(vc_analytical)-min(vc_analytical))/max(vc_analytical);
F2_vc1(z1)=sqrt(sum(f2_vc)/nx)/(mean(vc_analytical));

F1_hc1(z1)=(max(hc_analytical)-min(hc_analytical))/max(hc_analytical);
F2_hc1(z1)=sqrt(sum(f2_hc)/nx)/(mean(hc_analytical));

F1_dPmch1(z1)=(max(dPmch_analytical)-min(dPmch_analytical))/max(dPmch_analytical);
F2_dPmch1(z1)=sqrt(sum(f2_dPmch)/nx)/(mean(dPmch_analytical));

% fprintf('%3.2f %1.8f %1.3f %1.8f %1.8f %1.8f %3.1f %9.2f %9.2f %1.5f\n',Re_man2(i),Hchn(q,i)2(i),alp2(i),Wchn(q,i)2(i),Win2(i),Hman2(i),N2(i),h_BaseV2(i),P_pump_tot(i),F2_hc(i));
end
end
end
save('C:/Matlab Interface/signal1')




%
function sol=f_tur(x,q,k1,k2,k3,Rho,Mu,Dh_man,v_in,Beta1,Beta2,fRe_h)
% Part of the bvp4c solver (Governing equation)
% q(1)=v1 and q(2)=dv1/dx 
if q(1)<0
    q(1)=0;
end

if q(1)>1
    q(1)=1;
end

Re1=Rho*q(1)*Dh_man*v_in/Mu;
Re2=Rho*(1-q(1))*Dh_man*v_in/Mu;

if Re1>2300
fRe1=Re1*(0.79*log(Re1)-1.64)^-2/4;
else
    fRe1=fRe_h;
end

if Re2>2300
fRe2=Re2*(0.79*log(Re2)-1.64)^-2/4;
else
    fRe2=fRe_h;
end

v1=q(1);
ex1=q(2);
ex2=(-2*k1*q(2)+fRe2*k2*k1*(1-v1)-fRe1*k2*k1*v1-k1*(-Beta2+Beta2*q(1)-Beta1*q(1))*q(2))/(k3*q(2)-1);
sol=[ex1;ex2];


%
function sol=BC_tur(qw,qe,k1,k2,k3,Rho,Mu,Dh_man,v_in,Beta1,Beta2,fRe_h)
% Part of the bvp4c solver (The BC)
% qw=v(0)=1, qe=v(1)=0
sol=[qw(1)-1; qe(1)];
