classdef Cenario
	properties (Constant)
        mu = 398600.5;              % constante gravitacional da Terra
        rp=6371;                    % Raio do Planeta - km 
        rs = 695500;                % Raio do Sol - km
	end
     
	properties  ( Access = public )
        TempoTotalSimulacao;        % tempo total da simulação em segundos
        TempoIntervaloSimulacao;	% intervalo de simulação em segundos      
        TotalSatelites = 66;        % Total de Satélites da constelação
        TotalPlanos = 6;            % Total de Planos orbitais
        TotalSatPlanos = 11;        % Total de Satélites por plano
        BordaPolar = 75;            % Latitude da borda polar
        Zonas;                      % divisão da terra em zonas
        Constelacao;                % da constelação       
        Propagacao;                 % propagação no tempo
        Eclipse;                    % dados dos eclipses
    end
    
    methods
        
        %% Construtor
        function obj = Cenario(Data,TempoTotalSimulacao,TempoIntervaloSimulacao)
            % Constelação de satélites apartir do TLE
            obj.Constelacao	= obj.InfoTLE();
            % divisão da terra em zonas
            obj.Zonas =  Zonas();
            % tempo total da simulação em segundos
            obj.TempoTotalSimulacao = TempoTotalSimulacao;  
            % intervalo de simulação em segundos
            obj.TempoIntervaloSimulacao = TempoIntervaloSimulacao;	
            % Propagação dos satélites no Tempo
            obj.Propagacao = obj.Propagar(datevec(Data));  
            % Dados do eclipse
            obj.Eclipse = obj.DadosEclipse();
        end
        %% Extrai dados do ficheiro TLE
        function DadosOrbitais = InfoTLE(obj)        
                fd = fopen('TLE.txt', 'rb');
                A0 = fgetl(fd);
                A1 = fgetl(fd);
                A2 = fgetl(fd);
                DadosOrbitais = {};
                cont=1;
                while ischar(A2)          
                    % Epoca dos dados
                    DadosOrbitais(cont).epoca = obj.Epoca2Data(str2num(A1(19:32)));	 
                    % ecentricidade
                    DadosOrbitais(cont).e = str2double(A2(27:33))/ (1e7);
                    % inclinação
                    DadosOrbitais(cont).i = str2double(A2(9:16));  
                    % direção de nó ascendente
                    DadosOrbitais(cont).raan = str2double(A2(18:25)); 
                    % argumento de perigeu
                    DadosOrbitais(cont).w = str2double(A2(35:42));  
                    % anomalia Média
                    DadosOrbitais(cont).M = str2double(A2(44:51));    
                    % movimento médio
                    n = str2double(A2(53:63));  
                    DadosOrbitais(cont).n = n;   
                    % Semieixo maior
                    DadosOrbitais(cont).a = (obj.mu/(n*2*pi/(24*3600))^2)^(1/3);	
                    cont = cont + 1; 
                    A0 = fgetl(fd);
                    A1 = fgetl(fd);
                    A2 = fgetl(fd); 
                end
                fclose(fd);               
        end
        
        %% Retorna vetor de data/hora apartir da época do TLE
        function vetor_data = Epoca2Data(obj, tle_epoca )
            ymd=floor(tle_epoca);
            yr=fix(ymd/1000);
            dofyr=mod(ymd,1000);
            if (yr < 57)
                 year=  yr+ 2000;
            else
                 year=  yr+ 1900;
            end; 
            decidy=round((tle_epoca-ymd)*10^8)/10^8;
            temp=decidy*24;
            hh=fix(temp);
            temp=(temp-hh)*60;
            mm=fix(temp);
            temp=(temp-mm)*60;
            ss=floor(temp);
            nd = eomday(year,1:12);
            temp = cumsum(nd);
            month=find(temp>=dofyr, 1 );
            temp=temp(month)-dofyr;
            date=nd(month)-temp;
            vetor_data=[year,month,date,hh,mm,ss];
        end
 
        %% Retorna vetor de posição em ECI, apartir dos elementos orbitais
        function rECI = Oe2Eci(obj,SatNum,Data)  
                DeltaT = etime(Data,obj.Constelacao(SatNum).epoca);
                j2 = 1.0826359*10^-3;           % Pertuaão J2 (Achatamento da terra)             
                a = obj.Constelacao(SatNum).a;            % Semieixo maior       
                e = obj.Constelacao(SatNum).e;            % Ecentricidade       
                i = obj.Constelacao(SatNum).i;            % Inclinação    
                raan = obj.Constelacao(SatNum).raan;      % logitude do nodo ascendente      
                w = obj.Constelacao(SatNum).w;            % argumento de perigeu
                M = obj.Constelacao(SatNum).M;            % anomalia média          
                p = a*(1-e^2);                  % Parametro semi-latus rectum         
                T = 2*pi*sqrt(a^3/obj.mu);          % Periodo em segundos       
                Tp = (M/360*T);                 % Tempo de passagem pelo perigeu   
                ti = Tp+DeltaT;                 % Tempo Inicial da simulação considerando o perigeu    
                j2_raan = -(3/2)*(j2)*((obj.rp/p)^2)*sqrt(obj.mu/(a^3))*cosd(i);	   % Pertubação da longitude do nodo ascendente        
                j2_w = (3/4)*(j2)*((obj.rp/p)^2)*sqrt((obj.mu/a^3))*(5*cosd(i)^2-1) ;  % Pertubação do argumento de perigeu    
                M = mod(360*ti/T,360);          % Atualização da anomalia média  
                E = obj.CalcKeplerEq(M,e);      	% Atualização da anomalia excentrica              
                v = mod(atan2((sind(E)*(1-e^2)^.5),(cosd(E)-e)),2*pi)*180/pi;  % atualização da anomalia verdadeira       
                RAAN = ((raan*pi/180)+j2_raan*DeltaT)*180/pi;                  % atualização da longitude do nodo ascendete       
                w = ((w*pi/180)+j2_w*DeltaT)*180/pi;  % atualização do argumento de perigeu      
                Quaternio = quatnormalize(angle2quat(RAAN*pi/180,i*pi/180,v*pi/180+w*pi/180,'ZXZ'));
                % Calcula matriz de rotação do quaternio
                Rot_m = quat2rotm(Quaternio) ;  %   Matrix de rotação do quaternion 
                % raio orbital
                r = a*((1-(e*cosd(E))));
                % Calculo do vetor orbital 'r' em ECI
                rECI = Rot_m*[1;0;0]*r;    
        end   
        
        %% Retorna vetor de posição em ECEF
        function [rECEF] = Eci2Ecef(obj,rECI,Data)
            JD = juliandate(Data);
            Tj = (JD - 2451545.0)/36525.0;  
            GMST = mod(280.46061837 + 360.98564736629*(JD - 2451545) + 0.000387933*Tj^2- Tj^3/38710000,360);
            R3 = [cosd(GMST),sind(GMST),0;-sind(GMST),cosd(GMST),0;0,0,1];
            rECEF = R3*rECI; 
        end
        
        %%  Retorna Anomalia Ecêntrica - E
        function [E] = CalcKeplerEq(obj,M,e)
        tol = 10^-8;
        M = M*pi/180; % radianos
        Etemp = M;
        ratio = 1;
        while abs(ratio) > tol
            f_E = Etemp - e*sin(Etemp) - M;
            f_Eprime = 1 - e*cos(Etemp);
            ratio = f_E/f_Eprime;
            if abs(ratio) > tol
                Etemp = Etemp - ratio;
            else
                E = Etemp;
            end
        end
        E = mod(E,2*pi)*180/pi; % Graus
        end       
        
        %% calcula latitude e longitude a partir do vetor de posição em ECEF
        function [Lat,Long] = LatLong(obj,rECEF)
            r_delta = norm(rECEF(1:2));
            sinA = rECEF(2)/r_delta;
            cosA = rECEF(1)/r_delta;
            Long = atan2(sinA,cosA);
            if Long < -pi
                Long = Long + 2*pi;
            end
            Lat = asin(rECEF(3)/norm(rECEF));    
            Lat = Lat*180/pi;
            Long = Long*180/pi;
        end
        
        %% Calcula a distancia entre o centro de umbra e término de umbra e penumbra
        function [Eu Kp] = CalcEuKp(obj,vetor_rs,Sps)
           % Xu = Distância entre o Ápice do cone de umbra e o centro da terra - Km
            Xu = (obj.rp*2*Sps)/(obj.rs*2-obj.rp*2); 
            AlphaU = asin((2*obj.rp)/(2*Xu)); %ângulo de Umbra 
            Eu = (Xu-norm(vetor_rs))*tan(AlphaU);
            %Xp = Distância entre o Ápice do cone de peumbra e o centro da terra - Km
            Xp = (2*obj.rp*Sps)/(2*obj.rs+2*obj.rp);
            AlphaP = asin((2*obj.rp)/(2*Xp)); %ângulo de penumbra 
            Kp = (Xp+norm(vetor_rs))*tan(AlphaP);
        end 
        
        %% calcula ocorrência de eclipse
        function [EclipseStatus] = CalcEclipse(obj,Data,R)
            JD = juliandate(Data);
            %Julian centuries from 2000
            T = (JD - 2451545.0)/36525.0;  
            %geometric mean longitude of the sun
            L0 = mod(280.46645 + 36000.76983*T + 0.0003032*T^2,360);
            %mean anomaly of the sun
            M = mod(357.52910 + 35999.05030*T - 0.0001559*T^2 - 0.00000048*T^3,360);
            %eccentricity of the earth's orbit
            e = mod(0.016708617 - 0.000042037*T - 0.0000001236*T^2,360);
            %The equation of center for the sun, relative to its mean anomaly
            C = mod((1.914600-0.004817*T-0.000014*T^2)*sind(M)+(0.019993-0.000101*T)*sind(2*M)+0.000290*sind(3*M),360);
            %ecliptic longitude (lambda) of the sun
            Ls = mod((L0 + C),360)*pi/180;  
            %true anomaly of the sun
            f=mod(M+C,360);
            %earth-sun distance
            RS = ((1.000001018*(1 - e^2))/(1 + e*cosd(f)))*149597870.70;
            %Obliquity of the ecliptic
            Ep = mod(23 + 26/60 + 21.448/3600 - 46.8150/3600*T - (0.00059/3600)*T^2+ (0.001813/3600)*T^3,360)*pi/180;
            %Vetor UNITARIO sol   
            vetor_solar = [cos(Ls),sin(Ls)*cos(Ep),sin(Ls)*sin(Ep)];   
            vetor_rs = dot(R,vetor_solar)*vetor_solar;
            vetor_projecao = R'-vetor_rs; 
            [Eu Kp] = obj.CalcEuKp(vetor_rs,RS);
                % Eclipse status : 0 - sem eclipse, 1 = penumbra, 2 Umbra 
                if(dot(R,vetor_solar) < 0);   
                    NORM = norm(vetor_projecao);
                    if(NORM > Kp) 
                        EclipseStatus = 0;
                    end    
                    if(Eu<NORM && NORM <Kp)  
                        EclipseStatus = 1; % Penumbra                   
                    end   
                    if(NORM < Eu)            
                        EclipseStatus = 1; % Umbra
                    end 
                else
                    EclipseStatus = 0; 
                end   
        end      

        %% propagação dos satélites no tempo
        function Propagacao = Propagar(obj,Data,Interacoes)
            if nargin == 2
                Interacoes = ceil(obj.TempoTotalSimulacao/obj.TempoIntervaloSimulacao) ;
            end
            LoopId = 1;
            while LoopId <= Interacoes  
                fprintf('Cenário - Interação %d de %d\n',LoopId,Interacoes);  
                Propagacao(LoopId).Data = Data;
                Satelites = [];
               
                for SatNum=1:obj.TotalSatelites   
                    rECI  = obj.Oe2Eci(SatNum,Data);
                    Satelites(SatNum).rECI = rECI;
                    %## vetor de posição em ECEF
                    rECEF = obj.Eci2Ecef(rECI,Data);   
                    %##  determina a latitude e longitude
                    [Lat,Long] = obj.LatLong(rECEF);    
                    Satelites(SatNum).Latitude = Lat;          
                    Satelites(SatNum).Longitude = Long;
                    %## Ocorrência de Eclipse
                    EclipseStatus =  obj.CalcEclipse(Data,rECI);
                    % Situação do eclipse: 1 - em eclipse, 0 - Sol
                    Satelites(SatNum).EclipseStatus = EclipseStatus;     
                    % Tempo de Exposição em Eclipse
                    TempoEclipseIntervalo = 0;   
                    if EclipseStatus ~= 0 % Satélite inicia em eclipse
                        % verifica tempo que permanecerá em eclipse no intervalo
                        for i=1:obj.TempoIntervaloSimulacao 
                            %incremento de data
                            DataIntervalo = datevec(addtodate(datenum(Data), i, 'second'));
                            % calcula o raio orbital na nova data
                            R = obj.Oe2Eci(SatNum,DataIntervalo);
                            % calcula status
                            Eclipse = obj.CalcEclipse(DataIntervalo,R);  
                            if Eclipse ~= EclipseStatus % houve mudança
                                TempoEclipseIntervalo = TempoEclipseIntervalo+1;
                                break
                            else
                                TempoEclipseIntervalo = TempoEclipseIntervalo+1;
                            end
                        end
                    else % Satélite inicia no sol       
                        % verifica tempo que permanecerá no sol no intervalo
                        for i=1:obj.TempoIntervaloSimulacao 
                            %incremento de data
                            DataIntervalo = datevec(addtodate(datenum(Data), i, 'second'));
                            % calcula o raio orbital na nova data
                            R = obj.Oe2Eci(SatNum,DataIntervalo);
                            % calcula status
                            Eclipse = obj.CalcEclipse(DataIntervalo,R);  
                            if Eclipse ~= EclipseStatus % houve mudança                           
                                 TempoEclipseIntervalo = obj.TempoIntervaloSimulacao-i; 
                                 Satelites(SatNum).EclipseStatus = Eclipse; 
                                break
                            end
                        end
                    end
                     Satelites(SatNum).EclipseTempoIntervalo =  TempoEclipseIntervalo;   
                    if LoopId == 1 % condições iniciais 
                        Satelites(SatNum).EclipseTempoTotal	= TempoEclipseIntervalo;        
                    else
                        if(Satelites(SatNum).EclipseStatus ~=0) % Em Eclipse
                            %Acumula o tempo de eclipse da iteracao
                            %anterior com o tempo atual
                            Satelites(SatNum).EclipseTempoTotal = Propagacao(LoopId-1).Satelites(SatNum).EclipseTempoTotal+TempoEclipseIntervalo;
                            %fprintf(" SatNum %d Satelites(SatNum).EclipseTempoTotal %d \n", SatNum, Satelites(SatNum).EclipseTempoTotal);  
                        else
                            %Atualiza o tempo de eclipse para ser acumulado
                            Satelites(SatNum).EclipseTempoTotal = Propagacao(LoopId-1).Satelites(SatNum).EclipseTempoTotal; % Estava em zero(0);
                        end  
                     end     
                end
                Propagacao(LoopId).Satelites = Satelites;
                % Enlaces ativos
                q = obj.TotalSatPlanos; % numero de satelites por plano
                Adj = zeros(obj.TotalSatelites); % alocação matriz de adjacencia
                for k=1:obj.TotalPlanos  %varre os planos 
                    for j=1:q     %varre os satelites do plano            
                        SatNum = (k*q)-(q-j); %sequencia de numeração do satélite 1... t      
                        satAnt   = SatNum-1;  % satelite anterior da orbita
                        satPost  = SatNum+1;  % satelite posterior da orbita
                        planAnt  = SatNum-q;  % satelite do plano anteiror
                        planPost = SatNum+q;  % satelite do plano posterior 
                        %%%%% Adicionar Link para o GSO - Satpelite GeoEstacionário %%%%%% 
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if(j==1) satAnt   = k*q; end            % primeiro satelite anterior da orbita
                        if(j==q) satPost  = k*q-(q-1); end      % último satelite posterior da orbita       
                        if(k==1) planAnt  = ((obj.TotalPlanos-1)*q)+j; end    % primeiro satelite  da orbita anterior
                        if(k == obj.TotalPlanos) planPost = j; end              % ultimo satelite  da orbita anterior             
                        % ---- ISL do mesmo plano - Custo fixo            
                        La = sqrt(2)*norm(Satelites(SatNum).rECI)*sqrt(1-cosd(360/q)); % Custo do link 
                        Adj(SatNum,satAnt)  =  La; % Satelite anterior da mesma órbita
                        Adj(SatNum,satPost) =  La;  %Satelite posterior da mesma órbita          
                        if abs(Satelites(SatNum).Latitude) < obj.BordaPolar % se o satélite está em "Simulacao.BordaPolar", não possui ISL Intra-plano                   
                            Le = sqrt(2)*norm(Satelites(SatNum).rECI)*sqrt(1-cosd(360/(2*obj.TotalPlanos)))*cosd(Satelites(SatNum).Latitude); % Custo ISL inter plano 
                            if  abs(Satelites(planAnt).Latitude) < obj.BordaPolar % latitude do sat anterior                 
                               if(k>1) % seam 
                                     Adj(SatNum,planAnt) = Le;  % Satelite da órbita anterior
                               end
                            end            
                            if abs(Satelites(planPost).Latitude) < obj.BordaPolar % latitude do sat posterior
                                if(k<obj.TotalPlanos) %seam 
                                   Adj(SatNum,planPost) =  Le;  % Satelite da órbita posterior                                                                                
                                end
                            end 
                       end             
                    end
                end          
                
                Propagacao(LoopId).Enlaces = Adj;
                
                % SATÉLITES DAS RESPECTIVAS ZONAS E HORÁRIO LOCAL
                tam = size(obj.Zonas.LatLong);
                Sats = zeros(tam);
                Horas = zeros(tam);
                SatLatLong = [[Propagacao(LoopId).Satelites.Latitude]', [Propagacao(LoopId).Satelites.Longitude]'];
                diff = (360/tam(2))/2; % diferença de graus para considerar meridianos
                for i=1: tam(1)
                    for j=1:tam(2)           
                        % calcula o satélite mais próximo
                        latlon=obj.Zonas.LatLong(i,j);
                        [k, dist] = dsearchn([latlon{1}(1),latlon{1}(2)],SatLatLong);   
                        [dist, SatNum] = min(dist);
                        Sats(i,j) = SatNum;   
                        lon = latlon{1}(2); %longitude       
                        t = datenum(Data);
                        zd = timezone(lon+diff,'degrees');
                        t = addtodate(t, -zd, 'hour');
                        hora = hour(t);
                        if hora ==0 
                         hora = 24;
                        end
                        Horas(i,j) = hora;            
                    end
                end
                Propagacao(LoopId).Zonas.Satelites = Sats; % cobertura dos satélites por zona             
                % Fluxo da demanda de acordo com o perfil da hora 500TB/dia
                % No codigo orignal, não explica o valor de 500000000000000. Pesquisamos esse valor, que é 500TB por dia = 5x10^14 Bytes por dia.
                % O tamanho médio de um pacote enviado é de 210 bytes (Hussein, 2014).
                [ti,tj] = size(obj.Zonas.MatrizTrafego);
                Fk = zeros(max(obj.Zonas.Continentes(:))); % fluxo dos continentes
                Tij = zeros(ti,tj);
                for i=1:ti
                   for j=1:tj       
                       Tij(i,j) = (obj.Zonas.MatrizTrafego(i,j)/...
                           sum(obj.Zonas.MatrizTrafego(:)))*(500000000000000/3600)...
                           *((obj.Zonas.AtividadeHora(Horas(obj.Zonas.ZonaId(i))))/100);         
                       %continente de origem
                       CkOrg = obj.Zonas.Continentes(obj.Zonas.ZonaId(i));
                       %continente de destino
                       CkDst = obj.Zonas.Continentes(obj.Zonas.ZonaId(j));   
                       Fk(CkOrg,CkDst)  = Fk(CkOrg,CkDst)+Tij(i,j); 
                   end        
                end              
                %Fluxo entre continentes
                tk = max(obj.Zonas.Continentes(:));
                for k=1:tk
                    Propagacao(LoopId).Zonas.Fluxo(k,:) = Fk(k,:)/sum(Fk(k,:)); 
                end 
                % atualiza Data/hora da simulação em função do intervalo da simulação
                Data = datevec(addtodate(datenum(Data), obj.TempoIntervaloSimulacao, 'second'));                  
                % incrementa controle de loop
                LoopId = LoopId+1;       
            end
            %% Calcular Eclise total;
            
            
        end
        
        function Eclipse = DadosEclipse(obj)
            tam = size(obj.Propagacao,2); %% É o Número de Cenários (Adicionado por Renata Mota)
            Saida=zeros(obj.TotalSatelites,1);
            Entrada=zeros(obj.TotalSatelites,1);
            Tempo=zeros(obj.TotalSatelites,1);
            for s=1:obj.TotalSatelites
                   ciclo = 1; % Representa a Revolução (voltas ao redor da Terra) (Adicionado por Renata Mota)
                   for i=2:tam  
                       % fim do eclipse
                       % Acumula os tempos em eclipse  (Adicionado por Renata Mota)
                       if(obj.Propagacao(i).Satelites(s).EclipseStatus  == 0 && obj.Propagacao(i-1).Satelites(s).EclipseStatus ~=0)  
                           Saida(s,ciclo) = i-1;
                           EclipseTempo = obj.Propagacao(i-1).Satelites(s).EclipseTempoTotal ;  
                           Entrada(s,ciclo) = Saida(s,ciclo)-round(EclipseTempo/obj.TempoIntervaloSimulacao); 
                           Tempo(s,ciclo) = ceil(EclipseTempo/60);   
                           ciclo = ciclo+1;
                       end
                   end
            end
            Eclipse.Entrada = Entrada;
            Eclipse.Saida = Saida;
            Eclipse.Tempo = Tempo;
            Eclipse.TempoMaximo = max([Tempo]')';         
        end
    end    
end