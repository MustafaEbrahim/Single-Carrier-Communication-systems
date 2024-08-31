clc ; 
clear ; 
close all;
% parameters
QAM_symbol= 16 ;
No_Symbols = 1000 ;
NO_of_bits = QAM_symbol * No_Symbols * 3 ;
% noise
Noise_points = 19 ;
%Generate BER vectros
BPSK_BER_Actual  = zeros(Noise_points,1);
QPSK_BER_Actual  = zeros(Noise_points,1);
QPSK_BER_Actual_binary  = zeros(Noise_points,1);
Eight_BER_Actual = zeros(Noise_points,1);
QAM_BER_Actual   = zeros(Noise_points,1);
BFSK_BER_Actual = zeros(Noise_points,1);
BER_Actual = zeros(Noise_points,1);
    
BPSK_BER_theoritical = zeros(Noise_points,1);
QPSK_BER_theoritical = zeros(Noise_points,1);
Eight_PSK_BER_theoritical = zeros(Noise_points,1);
QAM_BER_theoritical = zeros(Noise_points,1);
BFSK_BER_theoritical = zeros(Noise_points,1);
    
No = zeros(1,Noise_points);
Eb_over_NO_dB= -4 : +1 : 14 ;
% test modulation and plot Constellation
Test_Modulation_Schemes();
% Run-all
for Schemes = 1 : 6 
    % Generating Data bits
    Data_bits_Generated = randi([0 1],1,NO_of_bits);  % generate the random logic bit-stream    
    if Schemes == 1
        disp("The current Modulation Scheme is BPSK ") ;
        bits_modulation = 1 ;
     elseif Schemes == 2
        disp("The current Modulation Scheme is QPSK ") ;
        bits_modulation = 2 ; 
     elseif Schemes == 3
        disp("The current Modulation Scheme is 8PSK ") ;
        bits_modulation = 3 ;
      elseif Schemes == 4
        disp("The current Modulation Scheme is 16-QAM ") ;
        bits_modulation = 4 ;
      elseif Schemes == 5
        bits_modulation = 1 ;
        disp("The current Modulation Scheme is BFSK ") ;
    else % QPSK Binary
        disp("The current Modulation Scheme is QPSK with Binary Encoding ") ;
        bits_modulation = 2 ;
    end 
    % Generate X_inphase and X_Quad vectros 
    X_inphase_Mapper = zeros (length (Data_bits_Generated)/bits_modulation ,1) ;
    X_Quad_Mapper    = zeros (length (Data_bits_Generated)/bits_modulation ,1) ;
   
    % Convert data bits into groups of bits
    dataBitsGrouped = reshape(Data_bits_Generated, bits_modulation , [])';
   % Calling Mapper Function 
    [X_inphase_Mapper , X_Quad_Mapper] = Mapper ( dataBitsGrouped  , Schemes );
    % Keep Symbol as X_i + j X_Q 
    Symbols = X_inphase_Mapper' + 1i.* X_Quad_Mapper';
    % Create constellation references 
    data_ref = 0 : pow2(bits_modulation)-1 ;
    % Calculate the binary representation
    data_test_binary = rem(floor(data_ref(:) * pow2(-bits_modulation+1:+1:0)), 2);
    [XI_ref , XQ_ref] = Mapper ( data_test_binary  , bits_modulation );  
    constellation_ref = XI_ref + 1i.* XQ_ref;
    
%% Channel
    % calculate E_avg 
    E_avg = sum(X_inphase_Mapper(:).*X_inphase_Mapper(:) + X_Quad_Mapper(:).*X_Quad_Mapper(:))/length (Symbols);
    for X= 1 : 19 
        Eb = 1 ; 
        No(X) = Eb/(10^(Eb_over_NO_dB(X)/10)); 
        % noise generated
        noise = randn(1 , length(Symbols))+1i .* randn(1 , length(Symbols)) ; 
        symbol_with_noise = Symbols + noise * (sqrt((No(X)/2) *(E_avg/bits_modulation))) ;
       if X == 15
            % Plotting the complex numbers as points without space in between
            figure;
            scatter(real(symbol_with_noise), imag(symbol_with_noise), 'filled','SizeData', 1);
            xlabel('Real Part');
            ylabel('Imaginary Part');
            title('constellation after adding noise due to AWGN channel ');
            grid on;
       end
        % call demapper function
        data_stream_demapper = Demapper( symbol_with_noise , constellation_ref , bits_modulation ,Schemes  );
        % check error 
        check_errors = (data_stream_demapper ~= Data_bits_Generated);
        BER_Actual(X) = length(check_errors(check_errors==1))/ length(data_stream_demapper) ;
        if Schemes == 1
            BPSK_BER_theoritical(X)= 0.5 *erfc (sqrt(1/No(X)));
            BPSK_BER_Actual(X) = BER_Actual(X);
        elseif Schemes == 2
            QPSK_BER_theoritical(X)= 0.5 *erfc (sqrt(1/No(X)));
            QPSK_BER_Actual(X) = BER_Actual(X);
        elseif Schemes == 3
            Eight_PSK_BER_theoritical(X)= (1/3) * erfc(sqrt(3*(1/No(X)))*sin(pi/8));
            Eight_BER_Actual(X) = BER_Actual(X);
        elseif Schemes == 4
            QAM_BER_theoritical(X)= (3/8)* erfc(sqrt(((1/No(X))/2.5)));
            QAM_BER_Actual(X) = BER_Actual(X);
        elseif Schemes == 5
            BFSK_BER_theoritical(X)= (0.5)* erfc(sqrt(0.5/No(X)));
            BFSK_BER_Actual(X) = BER_Actual(X);
        else % Schemes == 6
            QPSK_BER_theoritical(X)= 0.5 *erfc (sqrt(1/No(X)));
            QPSK_BER_Actual_binary(X) = BER_Actual(X);
        end 
    end
    
    figure ;
    semilogy (Eb_over_NO_dB,BER_Actual,'r','linewidth',1.5);
    hold on;
   if Schemes == 1
        semilogy(Eb_over_NO_dB,BPSK_BER_theoritical, '--','color',[0.5, 0.5, 0.5],'linewidth',2);
        title("BPSK BER Calculation");
        legend("Actual BER","Theoritical BER");
    elseif Schemes == 2
        semilogy(Eb_over_NO_dB,QPSK_BER_theoritical, '--','color',[0.5, 0.5, 0.5],'linewidth',2);
         title("Gray QPSK BER Calculation");
         legend("Actual BER","Theoritical BER");
    elseif Schemes == 3
        semilogy(Eb_over_NO_dB,Eight_PSK_BER_theoritical, '--','color',[0.5, 0.5, 0.5],'linewidth',2);
         title("8PSK BER Calculation");
         legend("Actual BER","Theoritical BER");
   elseif  Schemes == 4
        semilogy(Eb_over_NO_dB,QAM_BER_theoritical, '--','color',[0.5, 0.5, 0.5],'linewidth',2);
         title("16-QAM BER Calculation");
        legend("Actual BER","Theoritical BER");
   elseif  Schemes==5 
       semilogy(Eb_over_NO_dB,BFSK_BER_theoritical, '--','color',[0.5, 0.5, 0.5],'linewidth',2);
       title(" BFSK BER Calculation");
       legend("Actual BER","Theoritical BER");
   else 
        semilogy(Eb_over_NO_dB,QPSK_BER_theoritical, '--','color',[0.5, 0.5, 0.5],'linewidth',2);
         title("Binary and Gray QPSK  BER Calculation");
         hold on ;
         semilogy(Eb_over_NO_dB,QPSK_BER_Actual,'b','linewidth',1.5);
         title("QPSK Binary BER Calculation");
         legend("Binary QPSK BER","Theoritical BER " ,"Gray QPSK BER" );
    end 
    grid  on ;
    ylim([10^-4,10^1]);
    xlabel("Eb/N0 (db)");
    ylabel("BER");
    hold off;
end 
% plot all graph 
all_Simulated_BER = [BPSK_BER_Actual , QPSK_BER_Actual ...
                     QPSK_BER_Actual_binary , Eight_BER_Actual ,...
                     QAM_BER_Actual, BFSK_BER_Actual];
all_Theoritical_BER = [BPSK_BER_theoritical,Eight_PSK_BER_theoritical ,...
                         QAM_BER_theoritical ,BFSK_BER_theoritical ];              
figure;
colorOrder = get(gca, 'ColorOrder');
for C_S = 1 : 6
    semilogy(Eb_over_NO_dB,all_Simulated_BER(:,C_S)','Color', colorOrder(mod(C_S-1, size(colorOrder, 1)) + 1, :), 'LineWidth', 1.2);
    hold on;
    if C_S <= 4
        semilogy(Eb_over_NO_dB,all_Theoritical_BER(:,C_S)','--','Color', colorOrder(mod(C_S+5, size(colorOrder, 1)) + 1, :), 'LineWidth', 1.52);
        hold on;
    end 
end 
legend('Simulated BPSK  '       , 'theoritical BPSK ',...
       'Simulated QPSK (Gray)  ', 'theoritical 8PSK ',...
       'Simulated QPSK (Binary)', 'theoritical 16-QAM ',...
       'Simulated 8PSK '        , 'theoritical BFSK ',...
       'Simulated 16-QAM '      , 'Simulated BFSK ');
grid on
xlabel('Eb/No(dB)');
ylabel('Bit Error Rate ');
title('Bit Error Rate for all Schemes Modulation');
ylim([1e-4,1]);
xlim([-4 , 8]);

%% PSD of BFSK
% PSD Parameters
ensemble_size = 500;
num_bits_fsk = 100;
Tb = 5/1000  ; % 5 ms
num_samples_per_bit = 7 ;
sampling_time = Tb/num_samples_per_bit ; 
Total_samples = num_samples_per_bit * num_bits_fsk ; 
Eb_BFSK = 1 ;
step = (Tb/num_samples_per_bit);
delta_F = 1/(Tb) ;
% Generate time vector for One bit
t = (0:num_samples_per_bit-1)* step ;
%generating the baseband equivalent signal
S1_BB = sqrt(2*Eb_BFSK/Tb); %equivalent to 0  
S2_BB = sqrt(2*Eb_BFSK/Tb)*exp(1i*(2*pi*delta_F*t)); %equivalent to 1 
% Generate data matrix for X is No of Realization and y is num of bit +1 for delay
binary_sequence_PSD = randi([0,1], ensemble_size, num_bits_fsk+1); 
% Generate BFSK symbol
BFSK_symbols = zeros(ensemble_size,Total_samples+num_samples_per_bit);
for ensemble = 1 : ensemble_size
    for bit = 0 : num_bits_fsk
        if binary_sequence_PSD(ensemble , bit+1) == 1
            BFSK_symbols (ensemble,1+bit* num_samples_per_bit:(bit+1)*num_samples_per_bit)=S2_BB  ;
        else 
            BFSK_symbols (ensemble,1+bit* num_samples_per_bit:(bit+1)*num_samples_per_bit)=S1_BB ;
        end 
    end 
end 
% Generating Delay to add it
Td = randi([0,(num_samples_per_bit-1)],ensemble_size,1);
% Define matrix to store the values after adding the delay
Data_BFSK = zeros(ensemble_size, Total_samples);
% Apply the delay to the Binary BFSK Symbol
for i = 1:ensemble_size
    Data_BFSK(i,:)= BFSK_symbols(i, Td(i)+1 : Total_samples + Td(i));
end
% Define auto correlation matrix
initial_Stat_Auto_corr = zeros(size(Data_BFSK));
for i = 1:Total_samples
        for j = 1:Total_samples
            % Select two columns for element-wise multiplication 
            DOT_colmuns = conj(Data_BFSK(:, i)) .* Data_BFSK(:, j);
            % Perform element-wise multiplication
            if (j>=i)
                initial_Stat_Auto_corr(:,j-i+1) = initial_Stat_Auto_corr(:,j-i+1)+ DOT_colmuns ;
            end 
        end
end
Stat_Auto_corr = zeros ( 1 , Total_samples);
for i = 1 :Total_samples
    Stat_Auto_corr (1 , i ) = sum (initial_Stat_Auto_corr( : , i ) )/(Total_samples*ensemble_size) ;
end 
% Concatenate to get the final statistical autocorrelation
Stat_Auto_corr_full = cat (2, conj(fliplr(Stat_Auto_corr(2:end))), Stat_Auto_corr);
%Plotting PSD graph
delta = zeros(size(Stat_Auto_corr_full)); 
fs = 1 /sampling_time ;
M = -Total_samples+1 : Total_samples - 1;
f = M * fs * (Tb)/ (size(M,2))  ;  %Sampling frequency
delta(600) = 1/(4*Tb);  % Set delta function at desired frequency
delta(800) = 1/(4*Tb);
for H = 1 : length (f)
    delta(H) =abs( delta(H) + (4 * cos(pi * f(H))^2)/((pi)^2 *(4*(f(H))^2 -1)^2));
end 
PSD = abs(fftshift(fft(Stat_Auto_corr_full)))/(4*Total_samples); %normalized for 2Eb
figure('Name','PSD for BFSK Scheme');
plot(f , delta , 'r');
hold on ;
plot(f , PSD ,'b' );
legend('theoritical PSD ' ,'Simulated PSD  ' );
hold off;
title("PSD vs Frequency");
xlabel("Normalized frequency , F*T_b ");
ylabel("Normalized PSD ,S_B(F)/2E_b ");
ylim([0,0.8]);
xlim([-2,2]);
grid on;
%% Useful Functions 
% Mapper Block function
function [X_inphase_out , X_Quad_out] = Mapper ( Y_signal  , Schemes_used )
  X_inphase_out = zeros (length (Y_signal),1) ;
  X_Quad_out    = zeros (length (Y_signal),1) ;
  
  if Schemes_used == 1 
     [X_inphase_out , X_Quad_out] = BPSK_Mapper ( Y_signal ) ;
  elseif Schemes_used == 2 
     [X_inphase_out , X_Quad_out] = QPSK_Mapper ( Y_signal ) ;
  elseif Schemes_used == 3 
     [X_inphase_out , X_Quad_out] = Eight_PSK_Mapper ( Y_signal ) ;
  elseif  Schemes_used==4
     [X_inphase_out , X_Quad_out] = QAM_Mapper ( Y_signal ) ;
  elseif Schemes_used==5
     [X_inphase_out , X_Quad_out] = BFSK_Mapper ( Y_signal ) ;
  else %Schemes_used==6
     [X_inphase_out , X_Quad_out] = QPSK_Binary_Mapper ( Y_signal ) ;
  end 
end

% Demapper Block function
function data_Stream_Demapper = Demapper( Symbols_N , ref_symbol , No_Modulation_bits ,Schemes_used)
    XI_Correct = zeros (1,length (Symbols_N)) ;
    XQ_Correct = zeros (1,length (Symbols_N)) ;
    data_Demapper   = zeros (length (Symbols_N),No_Modulation_bits) ;
%     data_Stream_Demapper = zeros (1,length(Symbols_N)*No_Modulation_bits);
    % Get min distance
    for M = 1 : length (Symbols_N)
        [ value , min_Index ] = min (abs(Symbols_N(1,M) - ref_symbol ));
        XI_Correct(1,M) = real(ref_symbol(min_Index,1));
        XQ_Correct(1,M) = imag(ref_symbol(min_Index,1));
        % Convert to get the bits used
        if Schemes_used == 1 % BPSK 
            data_Demapper(M) = (XI_Correct(M)+1)/2 ;
         elseif Schemes_used==2
            data_Demapper(M,1) = (XI_Correct(1,M)+1)/2 ;
            data_Demapper(M,2) = (XQ_Correct(1,M)+1)/2 ;
        elseif Schemes_used == 3 % 8PSK
            phase = pi/4 ;
            i_Angle_used = angle(XI_Correct(M)+ 1i.*XQ_Correct(M)) / phase ;
            if i_Angle_used < 0
               i_Angle_used=  i_Angle_used + 8 ;
            end 
            dec_to_Binary = rem(floor(i_Angle_used * pow2(-No_Modulation_bits+1:+1:0)), 2);
            data_Demapper(M,:) = Gray_Generator (dec_to_Binary);
            
        elseif Schemes_used==4  % 16-QAM
            data_Demapper(M,1) = (sign(XI_Correct(M)) + 1) / 2 ;
            data_Demapper(M,2) = (3 - abs(XI_Correct(M)))  / 2 ;
            data_Demapper(M,3) = (sign(XQ_Correct(M)) + 1) / 2 ;
            data_Demapper(M,4) = (3 - abs(XQ_Correct(M)))  / 2 ;  
            
        elseif Schemes_used==5  % BFSK
            if real(Symbols_N(M)) > imag(Symbols_N(M))
                data_Demapper(M)= 0 ;
            else
                data_Demapper(M)= 1 ;
            end
        else   % Schemes_used==6   QPSK with binary encoding
            if [XI_Correct(1,M) XQ_Correct(1,M)] == [1 1]
                data_Demapper(M,:) = [1 0] ;
            elseif [XI_Correct(1,M) XQ_Correct(1,M)] == [1 -1]
                data_Demapper(M,:) = [1 1] ;
            elseif [XI_Correct(1,M) XQ_Correct(1,M)] == [-1 1]
                data_Demapper(M,:) = [0 1] ;
            else
                data_Demapper(M,:) = [0 0] ;
            end
        end
    end 
    % Convert reshaped matrix to a single column vector
     data_Stream_Demapper = reshape(data_Demapper',1, []);
end

% BPSK Mapping  
function [X_inphase , X_Quad] = BPSK_Mapper ( Y_signal_BPSK )
    X_inphase = zeros (length (Y_signal_BPSK),1) ;
    X_Quad    = zeros (length (Y_signal_BPSK),1) ; 
    X_inphase = 2*Y_signal_BPSK-1 ;
end

% QPSK Gray representation Mapping
function [X_inphase , X_Quad] = QPSK_Mapper ( Y_signal_QPSK )
    X_inphase = zeros (length (Y_signal_QPSK),1) ;
    X_Quad    = zeros (length (Y_signal_QPSK),1) ; 
     
    X_inphase = 2*Y_signal_QPSK(:,1)-1 ;
    X_Quad    = 2*Y_signal_QPSK(:,2)-1 ;
end

% 8PSK
function [X_inphase , X_Quad] = Eight_PSK_Mapper ( Y_signal_eight_PSK )
     X_inphase = zeros (length (Y_signal_eight_PSK),1) ;
     X_Quad    = zeros (length (Y_signal_eight_PSK),1) ;
     Z_complex = zeros (length (Y_signal_eight_PSK),1) ;
     % transfer form binary data to the binary of i_angle
     Binary_bits (:,1) = Y_signal_eight_PSK(:,1);
     for i = 2 : 3
         % XOR operation with shifted version of the sequence to get Binary
         Binary_bits(:,i) = xor ( Binary_bits(:,i-1) , Y_signal_eight_PSK(:,i) );
     end 
     i_Angle = Binary_bits * (pow2(2:-1:0)');
     
     for j = 1 : length (Z_complex)
         phase = pi/4 ;
         Z_complex (j) =  1* exp(i_Angle(j) * 1i * phase);
     end 
     X_inphase = real (Z_complex);
     X_Quad    = imag (Z_complex);
     
end 

% 16-QAM Mapping
function [X_inphase , X_Quad] = QAM_Mapper ( Y_signal_QAM )
 % [ Sign X_in , Value X_in ,Sign X_Quad , Value X_Quad ]
     X_inphase = zeros (length (Y_signal_QAM),1) ;
     X_Quad    = zeros (length (Y_signal_QAM),1) ;
     
     Mag_Inphase  = 3 - 2 * Y_signal_QAM(:,2) ;
     Sign_Inphase = 2 * Y_signal_QAM(:,1)-1 ; 
     
     Mag_Quad  = 3 - 2 * Y_signal_QAM(:,4) ;
     Sign_Quad = 2*Y_signal_QAM(:,3)-1 ; 
     
     X_inphase = Mag_Inphase .* Sign_Inphase  ;
     X_Quad =  Mag_Quad .* Sign_Quad ;
end

% BFSK Mapping
function [X_inphase , X_Quad] = BFSK_Mapper ( Y_signal_BFSK )
     X_inphase = zeros (length (Y_signal_BFSK),1) ;
     X_Quad    = zeros (length (Y_signal_BFSK),1) ;
     % when y_signal_BFSK = 0 ,mapper output = XI + j XQ = 1
     % while y_signal_BFSK = 1 ,  mapper output = XI + j XQ = j
     X_inphase =  1 - Y_signal_BFSK  ; % get the complemente of signal
     X_Quad = Y_signal_BFSK ;
end

% QPSK Binary Mapping
function [X_inphase , X_Quad] = QPSK_Binary_Mapper( Y_signal_QPSK_binary )
     X_inphase = zeros (length (Y_signal_QPSK_binary),1) ;
     X_Quad    = zeros (length (Y_signal_QPSK_binary),1) ;
 for k = 1 : length (Y_signal_QPSK_binary)
    if Y_signal_QPSK_binary(k,:) == [0 0]
        X_inphase(k) = -1 ;
        X_Quad(k)    = -1 ;
    elseif Y_signal_QPSK_binary(k,:)  == [0 1]
        X_inphase(k) = -1 ;
        X_Quad(k)    =  1 ;
        
    elseif Y_signal_QPSK_binary(k,:)  == [1 0]
        X_inphase(k) = 1 ;
        X_Quad(k)    = 1 ;
    else
        X_inphase(k) =  1 ;
        X_Quad(k)    = -1 ;
    end
 end 
end 

% Gray function Generator 
function Gray_bits = Gray_Generator ( Y_signal_Binary )
     Gray_bits = zeros (size (Y_signal_Binary)) ; 
     % Perform shift
     shiftedArray = [zeros(size (Y_signal_Binary,1),1), Y_signal_Binary(:,1:size (Y_signal_Binary,2) -1)];
     Gray_bits  = xor (Y_signal_Binary , shiftedArray) ;
end 

% test Modulations
function Test_Modulation_Schemes()
disp("Test Modulation Schemes is running... ") ;
 for k = 1 : 6
    if k==1 % BPSK
      num_bits = 1;
    elseif k ==2 %Gray QPSK
      num_bits = 2;
    elseif k==3 % 8PSK
      num_bits = 3;
    elseif k==4 % 16-QAM
      num_bits = 4;  
    elseif k==5 % BFSK
      num_bits = 1;
    else % K == 6 Binary QPSK
      num_bits = 2;  
    end 
    data_test = 0 : pow2(num_bits)-1 ;
    % Calculate the binary representation
    data_test_binary = rem(floor(data_test(:) * pow2(-num_bits+1:+1:0)), 2);
    % rem () take reminder of divsion
    % if num_bits = 4
    % data_test(:) convert from Row to column vectors
    % 15 * (2^-3 , 2^-2 , 2^-1 , 2^0 ) = ( 1 ,3 ,7 , 15 )
    % rem (( 1 ,3 ,7 , 15 ) / 2) = ( 1 ,1 ,1 ,1 ) 
    [X_inphase_Mapper , X_Quad_Mapper] = Mapper ( data_test_binary  , k );
    Symbols_constellation = X_inphase_Mapper + 1i.* X_Quad_Mapper;
    % Plot the complex number
    figure;
    hold on;
    for j = 1:length(Symbols_constellation)
        % Plot complex number
        plot(real(Symbols_constellation(j)), imag(Symbols_constellation(j)),'bo' );
        % Add text annotation for magnitude
        text(real(Symbols_constellation(j)), imag(Symbols_constellation(j)),...
        [' ' num2str(data_test_binary(j,1:num_bits))], 'VerticalAlignment',...
        'bottom' , 'FontSize', 10);
    end
    hold on;
    xlabel('Real Part');
    ylabel('Imaginary Part');
    ylim([-2,2]);
    xlim([-2,2]);
    grid on;
    if k == 1
        title('Constellation of BPSK Modulation ');
    elseif k==2
        title('Constellation of Gray QPSK Modulation ');
        % Define square vertices
        x = [1  1 -1 -1  1 ];
        y = [1 -1 -1  1  1 ];
        % Plot square
        plot(x, y, 'k', 'LineWidth', 0.5);
        axis equal; 
    elseif k == 3
        title('Constellation of 8PSK Modulation ');
        % Define angles
        theta = linspace(0, 2*pi, 100); % 100 points around the circle
        % Define radius
        r = 1;
        % Calculate circle coordinates
        x = r * cos(theta);
        y = r * sin(theta);
        % Plot circle
        plot(x, y,'k','LineWidth', 0.5);
        axis equal; 
    elseif k == 4
        title('Constellation of 16-QAM Modulation ');
        ylim([-4,4]);
        xlim([-4,4]);
    elseif k == 5
        title('Constellation of BFSK Modulation ');
    else
       title('Constellation of Binary QPSK Modulation ');
        % Define square vertices
        x = [1  1 -1 -1  1 ];
        y = [1 -1 -1  1  1 ];
        % Plot square
        plot(x, y, 'k', 'LineWidth', 0.5);
        axis equal;
    end 
    hold off;
 end 
 disp("Testing finshed... ") ;
end 
