classdef Simular

     properties
        Cenario;            % Dados do cenário da simulação
        TotalFontes;        % números de terminais para simulação
        CBRs;               % taxas de bits para simulação
        Metricas;           % Métricas da simulação
        PotenciaCP;         % consumo de potência operação nominal
        PotenciaTX;         % consumo de potência na transmissão
        PotenciaRX;         % consumo de potência na recepção
        PotenciaON;         % consumo de potência na operação normal
        PotenciaCG;         % potência de captação de energia
        ISLCap;             % Capacidade dos links   
        Resultados;         % variável de resultados
        Lb;
        Lc;
        TipoAdaptacao;
        Pesos;
        PesosLimiar;
        constanteA;
    end
    
    methods
        
        % construtor
        function obj = Simular(Data,TempoTotalSimulacao,TempoIntervaloSimulacao)
          % Gera cenário da simulação
          obj.Cenario = Cenario(Data,TempoTotalSimulacao,TempoIntervaloSimulacao);
        end

        function obj = set.Cenario(obj,val)
          obj.Cenario = val;
        end
        
        function obj = set.TotalFontes(obj,val)
          obj.TotalFontes = val;
        end     
        
        function obj = set.CBRs(obj,val)
          obj.CBRs = val;
        end  
        
        function obj = set.Metricas(obj,val)
          obj.Metricas = val;
        end  
        
        function obj = set.PotenciaCP(obj,val)
          obj.PotenciaCP = val;
        end   
        
        function obj = set.PotenciaON(obj,val)
          obj.PotenciaON = val;
        end 
        
        function obj = set.PotenciaTX(obj,val)
          obj.PotenciaTX = val;
        end        
        
        function obj = set.PotenciaRX(obj,val)
          obj.PotenciaRX = val;
        end
        
        function obj = set.PotenciaCG(obj,val)
          obj.PotenciaCG = val;
        end  
        
        function obj = set.ISLCap(obj,val)
          obj.ISLCap = val;
        end        
        
        function obj = set.Resultados(obj,val)
          obj.Resultados = val;
        end
        
        function obj = set.Lc(obj,val)
          obj.Lc = val;
        end
        
        function obj = set.Lb(obj,val)
          obj.Lb = val;
        end
        
        function Resultados = Roteamento(obj,loopId,Interacoes)
                  
            if nargin == 2
                Interacoes = size(obj.Cenario.Propagacao,2);
            end   
            Resultados = {};
            Edges = {};
            % ---------------- simulação no tempo -------------------------------------
            LoopId = 1;
            
            %fprintf(' Teste -- Simulação %d (Roteamento %d/%d)\n',loopId,LoopId,Interacoes);  
            
            while LoopId <= Interacoes 
                fprintf('Simulação %d (Roteamento %d/%d)\n',loopId,LoopId,Interacoes);       
                % ---------------------- Geração Fontes / Destino -----------------   
                % Acesso: Resultados(TerminalID).Fontes 
                % pares - seleção dos continentes de origem
                %CkSrc = randsrc(obj.TotalFontes,1,[1:max(obj.Cenario.Zonas.Continentes(:)); ...
                %obj.Cenario.Zonas.ContinentesHosts./sum(obj.Cenario.Zonas.ContinentesHosts(:))]);
                CkSrc = randsrc(obj.TotalFontes,1,1:max(obj.Cenario.Zonas.Continentes(:)));
                for fonte_id=1: obj.TotalFontes
                    % continente de origem
                    Pares(fonte_id).ContinenteOrigem = CkSrc(fonte_id); 
                    % id das zonas de origem com demanda de tráfego
                    IdCkSrc = find(obj.Cenario.Zonas.Continentes...
                        .*(obj.Cenario.Zonas.DensidadeHosts>0) == CkSrc(fonte_id));
                    % seleção da zona de origem
                    %ZkSrc = randsrc(1,1,[IdCkSrc'; obj.Cenario.Zonas.DensidadeHosts(IdCkSrc)'...
                     %   ./sum(obj.Cenario.Zonas.DensidadeHosts(IdCkSrc))]);  
                    ZkSrc = randsrc(1,1,IdCkSrc');  
                    
                    Pares(fonte_id).ZonaOrigem = ZkSrc; 

                    % Satélites do continente de origem
                    SatSrc = obj.Cenario.Propagacao(LoopId).Zonas.Satelites(ZkSrc)';
                    Pares(fonte_id).SateliteOrigem = SatSrc;
                    % Continente de destino
                    CkDst = randsrc(1,1,[1:max(obj.Cenario.Zonas.Continentes(:));...
                        obj.Cenario.Propagacao(LoopId).Zonas.Fluxo(CkSrc(fonte_id),:)]);                  
                    Pares(fonte_id).ContinenteDestino = CkDst;      
                    % seleção da zona/satélite de destino - evita selecionar mesmo
                    while(1)
                        IdCkDst = find(obj.Cenario.Zonas.Continentes...
                            .*(obj.Cenario.Zonas.DensidadeHosts>0) == CkDst);
                        %ZkDst = randsrc(1,1,[IdCkDst'; obj.Cenario.Zonas.DensidadeHosts(IdCkDst)'...
                         %   ./sum(obj.Cenario.Zonas.DensidadeHosts(IdCkDst))]);  
                        ZkDst = randsrc(1,1,IdCkDst'); 
                        Pares(fonte_id).ZonaDestino = ZkSrc;
                        SatDst = obj.Cenario.Propagacao(LoopId).Zonas.Satelites(ZkDst)';
                        Pares(fonte_id).SateliteDestino = SatDst;
                        % evita selecionar mesmos atélite para origem e destino
                        if SatSrc ~= SatDst
                            break;
                        end       
                    end
                end
                Resultados(LoopId).Fontes = Pares; 
                % ----------------------- Enlaces Ativos ------------------------------
                % Matrix com o tempo de propagação de cada enlace - milisegundos
                MatDelay = ((obj.Cenario.Propagacao(LoopId).Enlaces...
                ./299792458).*1000000);   %%  velocidade da luz no vácuo = 299.792.458 m/s
                %% Tempo = Enlace (km)/ V (m/s) = Enlace*10^3 m / 299.792.458 m/ 10^3 milissegundos 
                %% Tempo = Enlace (km)/ V (m/s) = Enlace*10^3 * 10^3 milissegundos / 299.792.458
                %% Tempo = Enlace (km)/ V (m/s) = Enlace*1000000 milissegundos / 299.792.458
               
                % Matriz de adjacencia da rede
                MatRede  = (MatDelay>0); % Matriz Binaria com  ={ 1, para  MatDelay>0; 0, caso contrário} 
          
                % --------------------- Taxas de Bits --------------------------------
                for cbr_id = 1: size(obj.CBRs,2)
                    cbr = obj.CBRs(cbr_id);
                    % --------------------- Métricas ----------------------------------
                    for metrica_id=1:size(obj.Metricas,1)
                        Metrica = num2str(cell2mat(obj.Metricas(metrica_id)));   
                        % --------------- matriz de energia ---------------------------
                         ConsumoEnergia = zeros(obj.Cenario.TotalSatelites,1);
                         CicloVida = zeros(obj.Cenario.TotalSatelites,1);
                         GD = zeros(obj.Cenario.TotalSatelites,1);
                         DOD = zeros(obj.Cenario.TotalSatelites,1);
                         DOD_ant = zeros(obj.Cenario.TotalSatelites,1);
                         
                        if LoopId == 1 % ######### Condições iniciais  
                            EnergiaInicio = ones(obj.Cenario.TotalSatelites,1)*obj.PotenciaCP;
                            DOD = zeros(obj.Cenario.TotalSatelites,1);
                            %ConsumoEnergia = zeros(obj.Cenario.TotalSatelites,1);
                        else % ######### Condições não iniciais
                            % Energia inicio do intervalo de tempo é igual a 
                            % energia residual do final do intervalo anterior  
                            EnergiaInicio = Resultados(LoopId-1).(Metrica).EnergiaFinal(:,cbr_id);
                            DOD_ant =  Resultados(LoopId-1).(Metrica).DOD(:,cbr_id); %DOD da Iteracao anterior
                            %ConsumoEnergia = Resultados(LoopId-1).(Metrica).ConsumoEnergia(:,cbr_id);
                        end 
                        Resultados(LoopId).(Metrica).EnergiaInicio(:,cbr_id) = EnergiaInicio;
                        
                        Resultados(LoopId).(Metrica).ConsumoEnergia(:,cbr_id) = ConsumoEnergia;
                        
                        % Matriz de capacidade de Energia já desconta a energia de 
                        % Operação nominal dos satélites
                        EnergiaON = ones(obj.Cenario.TotalSatelites,1)...
                            .*obj.PotenciaON.*obj.Cenario.TempoIntervaloSimulacao;
                        %%% Energia residual corrente dos satelites %%%%%
                        Energia =  EnergiaInicio  -  EnergiaON;
                        % Atualiza Consumo de Energia
                        ConsumoEnergia = EnergiaON; 
                        
                        % --- Matrizes de Capacidade, Demanda e Demanda Atendida ------
                        % Capacidade dos enlaces
                        Capacidade = MatRede.*obj.ISLCap; % em Mbps
                        % Dados trafegados por enlace
                        MatTrafego = MatRede.*0;
                        % Demanda Origem/Detino
                        MatDemanda = MatRede.*0;       
                        % DemandaAtendida Origem/Detino
                        MatDemandaAtendida = MatRede.*0;
                       
                        % ----  Varre os pares de terminais (origem e destino) --------  
                        for fonte_id=1: obj.TotalFontes   
                            Origem  = Resultados(LoopId).Fontes(fonte_id).SateliteOrigem;
                            Destino = Resultados(LoopId).Fontes(fonte_id).SateliteDestino;
                            Demanda = cbr; % CBR em Mbps
                            % --------------------- Aplicação de Métricas -------------
                            switch Metrica    
                                case 'TP'    
                                    G = digraph(MatDelay);
                                    % tempo propagação                     
                                    Tmin = min(G.Edges.Weight); % utopia   
                                    Tmax = max(G.Edges.Weight); % nadir   
                                    % tempo de propagação normalizado
                                    TPNorm = (G.Edges.Weight - Tmin)/(Tmax - Tmin);
                                    G.Edges.Weight = TPNorm;           
                                case 'LASER'        
                                    G = digraph(MatDelay);
                                    % tempo propagação                     
                                    Tmin = min(G.Edges.Weight); % utopia   
                                    Tmax = max(G.Edges.Weight); % nadir   
                                    % tempo de propagação normalizado
                                    TPNorm = (G.Edges.Weight - Tmin)/(Tmax - Tmin);
                                    [Si,Sj] = findedge(G);
                                    Ei = [obj.Cenario.Propagacao(LoopId).Satelites(Si).EclipseStatus]';
                                    Ej = [obj.Cenario.Propagacao(LoopId).Satelites(Sj).EclipseStatus]'; 
                                    %Normalização do Nível de bateria;
                                    Bi = Energia(Si)./obj.PotenciaCP;
                                    Bj = Energia(Sj)./obj.PotenciaCP;
                                    Dij = Ei./Bi + Ej./Bj;
                                    Dij(isnan(Dij)) = 0; 
                                    Dmin = min(Dij); % UTOPIA   
                                    Dmax = max(Dij); % NADIR  
                                    DijNorm = (Dij-Dmin)/(Dmax-Dmin);
                                    G.Edges.Weight = 0.5.*TPNorm+0.5.*DijNorm;   
                                    
                                case 'PROPOSTA' 
                                    G = digraph(MatDelay);
                                    % tempo propagação                     
                                    Tmin = min(G.Edges.Weight); % utopia   
                                    Tmax = max(G.Edges.Weight); % nadir     
                                    % tempo de propagação normalizado
                                    TPNorm = (G.Edges.Weight - Tmin)/(Tmax - Tmin);
                                    [Si,Sj] = findedge(G);
                                    Ei = [obj.Cenario.Propagacao(LoopId).Satelites(Si).EclipseStatus]';
                                    Ej = [obj.Cenario.Propagacao(LoopId).Satelites(Sj).EclipseStatus]'; 
                                    Bi = Energia(Si)./obj.PotenciaCP;
                                    Bj = Energia(Sj)./obj.PotenciaCP;
                                    % Nova modelagem, considerando apenas o
                                    Peso = ones(size(G.Edges.Weight,1),3);
                                    Peso(:,1) = obj.Pesos(1);
                                    Peso(:,2) = obj.Pesos(2);
                                    Peso(:,3) = obj.Pesos(3);
                                    
                                    xij = zeros(size(G.Edges.Weight,1),1);
                                    pos = zeros(size(G.Edges.Weight,1),1);
                                    if ( obj.Lb > 0)
                                        BijBmax = [Bi(:,:) Bj(:,:)];
                                        xij = any (BijBmax < obj.Lb, 2);
                                        xij = double(xij(:));
                                        pos = find( xij == 1); % posicao dos satelites a serem removidos
                                    end
                                      
                                    if (obj.TipoAdaptacao == 2)
                                       Dij =  -(Bi + Bj); 
                                    else
                                       Dij = Ei./Bi + Ej./Bj; 
                                    end
                                    
                                    Cij = zeros(size(G.Edges.Weight,1),1);
                                    CijResidual = zeros(size(G.Edges.Weight,1),1);
                                    for s_ij = 1:size(Cij,1)
                                        CijResidual(s_ij,1) = Capacidade(Si(s_ij),Sj(s_ij)) / obj.ISLCap ;
                                        Cij(s_ij,1) = 1 - CijResidual(s_ij,1); 
                                        %Cij(s_ij,1) = 1-(Capacidade(Si(s_ij),Sj(s_ij))/obj.ISLCap);
                                    end 
                                    
                                    % Remove dos enlaces - Os satélites e
                                    % links da bateria que não atendem a
                                    % capacidade da bateria e capacidade do link
                                    
                                    Zij = zeros(size(G.Edges.Weight,1),1);
                                    posZij = zeros(size(G.Edges.Weight,1),1);
                                    if ( obj.Lc > 0)
                                      Zij = any (CijResidual(:,:) < obj.Lc, 2); 
                                      Zij = double(Zij(:));
                                      posZij = find(  Zij == 1);
%                                       pos = union(pos, posZij)
                                    end
                                    
                                    if (sum (pos) > 0 && sum (posZij) > 0 )
                                         pos = union(pos, posZij);
                                    elseif ( sum (pos) == 0 && sum (posZij) > 0)
                                        pos = posZij;
                                    end
                                    
                                    if (sum (pos) > 0 )
                                         Cij(pos, :) = NaN; % Atribui NaN para remover da normalização;
                                         Dij(pos, :) = NaN; 
                                         Peso(pos,1) = obj.PesosLimiar(1);
                                         Peso(pos,2) = obj.PesosLimiar(2);
                                         Peso(pos,3) = obj.PesosLimiar(3);
                                    else
                                         Dij(isnan(Dij)) = 0;
                                    end
                                    
                                    Dmin = min(min(Dij,[],'omitnan')); % utopia
                                    Dmax = max(max(Dij,[],'omitnan')); % nadir 
                                    DijNorm = (Dij-Dmin)/(Dmax-Dmin);

                                    Cmin = 0;
                                    Cmax = 1;
                                    CijNorm = (Cij - Cmin)./(Cmax-Cmin);
                                    
                                    if (sum (pos) > 0 )
                                      CijNorm(pos, :) = inf; % Remove os links sem capacidade;
                                      DijNorm(pos, :) = inf; % Remove Arestas com Satélite Sem bateria;  
                                      TPNorm(pos, :) = inf;  % Remove Arestas com Satélite Sem bateria;  
                                      %Se for NaN, Not a Number, atribui
                                      %infinito..
%                                       CijNorm(isnan(CijNorm)) = inf; 
%                                       DijNorm(isnan(DijNorm)) = inf ; 
%                                       TPNorm(isnan(TPNorm)) = inf; 
                                    end
                                      % Se não for um número desabilida a
                                      % rota;
                                     CijNorm(isnan(CijNorm)) = 0; 
                                     DijNorm(isnan(DijNorm)) = 0 ; 
                                     TPNorm(isnan(TPNorm)) = 0;  
                                    
                                    %%%%METRICA PROPOSTA : w1, w2, w3 %%%%%%%
                                    %G.Edges.Weight = 0.35.*TPNorm+0.35.*DijNorm+0.3.*CijNorm;
                                    TDC = [TPNorm DijNorm CijNorm];
                                    G.Edges.Weight = diag( Peso(:,:)  * TDC(:,:)');
                                    %Resultados(LoopId).(Metrica).Fontes.demanda(fonte_id,cbr_id) = Demanda;
                             end
                            
                             DemandaAtendida = Demanda; % em Mbps - Megabits por segundo
                             % matrix de demandas origem/destino no intervalo do tempo
                              MatDemanda(Origem,Destino) = MatDemanda(Origem,Destino) ...
                                + Demanda*obj.Cenario.TempoIntervaloSimulacao;  
                              
                              [rota, custo] = shortestpath(G,Origem,Destino); %%%%% Esse Custo Não é Utilizado %%%%%%
                                % Energia dos satélites da rota
                                EnergiaRota =  Energia(rota);  
                                % Energia  mínima da rota
                                
                                MinEnergiaRota = min(EnergiaRota);
                                %%%% Energia Tx, Rx é por enlace entre os satélites k-l
                                % Energia Requerida para transmissão
                                EnergiaTx = ((DemandaAtendida/obj.ISLCap)*obj.PotenciaTX)*obj.Cenario.TempoIntervaloSimulacao;
                                % Energia Requerida para recepcao
                                EnergiaRx = ((DemandaAtendida/obj.ISLCap)*obj.PotenciaRX)*obj.Cenario.TempoIntervaloSimulacao;
                                %Total de energia requerida para transmissão e recepção
                                EnergiaTxRx = EnergiaTx + EnergiaRx;

                                % verifica se há energia suficente para transmitir
                                if MinEnergiaRota < EnergiaTxRx
                                   DemandaAtendida = 0;
                                end  
                                
                            CapRota = [];
                            % total de saltos da rota
                            saltos = size(rota,2)-1;
                            %Encontrou uma rota com capacidade;
                            for st=1:saltos        
                                % capacidade do enlace do salto
                                CapEnlace = Capacidade(rota(st),rota(st+1)); % Mbps
                                CapRota = cat(1,CapRota,CapEnlace);
                                if CapEnlace < DemandaAtendida % Se ainda tem capacidade do Enlace
                                   % limita a demanda à capacidade
                                   DemandaAtendida = CapEnlace;
                                end
                            end 
                            
                            % tempo de propagação do enlace em milisegundos
                            Delay = 0;
                            % Caso haja dados para trafegar na rota efetua o roteamento
                            if DemandaAtendida > 0
                                
                                for st=1:saltos        
                                    % atualiza energia - transmissão
                                    Energia(rota(st)) = Energia(rota(st)) - EnergiaTx;
                                    % atualiza energia - recepção
                                    Energia(rota(st+1)) = Energia(rota(st+1)) - EnergiaRx; 
                                    % atualiza a Capacidade do enlace
                                    Capacidade(rota(st),rota(st+1))= Capacidade(rota(st),rota(st+1))...
                                        - DemandaAtendida;
                                    % Matriz de tráfego dos enlaces
                                    MatTrafego(rota(st),rota(st+1))...
                                        = MatTrafego(rota(st),rota(st+1))+DemandaAtendida...
                                        *obj.Cenario.TempoIntervaloSimulacao; %% Tráfego em Megabits(Mb) = Demanda(MB/s)* TempoIntervaloSimulacao (s);
                                    % Delay da rota (tempo de propagação)
                                    Delay = Delay+MatDelay(rota(st),rota(st+1));   
                                    
                                    %%%%%####Soma a energia Transmissão  e
                                    %%%%%Recepcao Gastas;
                                    ConsumoEnergia(rota(st)) = ConsumoEnergia(rota(st)) + EnergiaTx; % Consumo de Energia Transferência
                                    ConsumoEnergia(rota(st+1)) =  ConsumoEnergia(rota(st+1)) + EnergiaRx; % Consumo de Energia de Recepção
                                    %%%%%%%
                                end
                            else
                                % Não houve entrega dos dados
                                saltos = 0;          
                            end
                             % Matriz de demanda atendida no intervalo de tempo
                             MatDemandaAtendida(Origem,Destino) = MatDemandaAtendida(Origem,Destino)...
                                + DemandaAtendida*obj.Cenario.TempoIntervaloSimulacao; 
                             % Resumo Fontes
                             Resultados(LoopId).(Metrica).Fontes.rota(fonte_id,cbr_id) = {rota};
                             Resultados(LoopId).(Metrica).Fontes.demanda(fonte_id,cbr_id) = Demanda;
                             Resultados(LoopId).(Metrica).Fontes.capacidadesrota(fonte_id,cbr_id) = {CapRota};
                             Resultados(LoopId).(Metrica).Fontes.enlacessaturados(fonte_id,cbr_id) = size(find(CapRota==0),1);
                             Resultados(LoopId).(Metrica).Fontes.demandaatendida(fonte_id,cbr_id)...
                                 = DemandaAtendida;
                             Resultados(LoopId).(Metrica).Fontes.saltos(fonte_id,cbr_id) = saltos;
                             %if Delay ==0
                                % Delay = inf;
                             %end
                             Resultados(LoopId).(Metrica).Fontes.tempopropagacao(fonte_id,cbr_id) = Delay;   
                             Resultados(LoopId).(Metrica).Fontes.eclipserota(fonte_id,cbr_id)...
                                 = size(find([obj.Cenario.Propagacao(LoopId).Satelites(rota).EclipseStatus] ==1),2); 
                            
                            %%%%%Consumo de energia por enlace do satélite k %%%% 
%                           Consumo =  Resultados(LoopId).(Metrica).ConsumoEnergia(Origem,cbr_id);
%                           Consumo =  EnergiaON(Origem) + TotalTxRx; %+ Consumo;
%                           Resultados(LoopId).(Metrica).ConsumoEnergia(Origem,cbr_id) = Consumo;                             
                        end
                      
                        % demandas
                        %Resultados(LoopId).(Metrica).Demanda(:,cbr_id) = {MatDemanda};
                        %Resultados(LoopId).(Metrica).DemandaAtendida(:,cbr_id) = {MatDemandaAtendida};     
                        % Capacidade dos enlaces ao final do intervalo
                        Resultados(LoopId).(Metrica).Trafego(:,cbr_id) = {MatTrafego};  % Em Megabytes;           
                        % Capacidade dos enlaces ao final do intervalo
                        Resultados(LoopId).(Metrica).Capacidade(:,cbr_id) = {Capacidade};
                        % Captação de energia no intervalo
                        CaptacaoEnergetica = obj.PotenciaCG*max(0,obj.Cenario.TempoIntervaloSimulacao-[obj.Cenario.Propagacao(LoopId).Satelites.EclipseTempoIntervalo]');           
                        %%%%%Captação Energética ...
                        %%% Adicionada para Cálculo da Captação Energética
                        %CaptacaoEnergetica = Resultados(LoopId).(Metrica).CaptacaoEnergetica(:,cbr_id) + CaptacaoEnergetica
                        %%%%Calculo origninal %%%
                        Resultados(LoopId).(Metrica).CaptacaoEnergetica(:,cbr_id) = CaptacaoEnergetica;
                        % Atualiza matriz energética
                        
                        %%Energia Residual = Energia 
                          A = obj.constanteA;
                        % Energia(:) =  Energia(:) - Energia(:)*0.25; % teste diminui 25%
                          Dt2 =  (obj.PotenciaCP - Energia(:) ) ./ obj.PotenciaCP;
                          DOD(:,:) =  Dt2(:,:);
                          Dt1 = DOD_ant(:,:);
                         
                         for k=1:obj.Cenario.TotalSatelites
                            if( Dt2(k) > Dt1(k))
                                CicloVida(k,:) = Dt2(k) * 10^(A*(Dt2(k)-1)) - Dt1(k) * 10^(A*(Dt1(k)-1));
                            else
                                CicloVida(k,:) = 0;
                            end
                            GD(k)= Dt2(k) * 10^(A*(Dt2(k)-1));
                         end
                        
                        Energia = min(obj.PotenciaCP, Energia + CaptacaoEnergetica);
                        Energia = max(0, Energia); % Adicionado para evitar Energia com valor negativo (por Renata Mota)
                       
                        Resultados(LoopId).(Metrica).EnergiaFinal(:,cbr_id)  = Energia;
                        Resultados(LoopId).(Metrica).ConsumoEnergia(:,cbr_id) = ConsumoEnergia;
                        
                        Resultados(LoopId).(Metrica).DOD(:,cbr_id) = DOD;
                        Resultados(LoopId).(Metrica).CicloVida(:,cbr_id) = CicloVida;
                        Resultados(LoopId).(Metrica).GD(:,cbr_id)  = GD;
                       
                    end 
                end   
                LoopId = LoopId+1;
            end                       
        end   
    end   
end