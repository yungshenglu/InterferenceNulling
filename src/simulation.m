clear;
close all;
clc;

% Define simulated SNR value
SNR_simulated = [ 5 , 10 , 15 , 20 , 25 ];
SNR_decoded = [];


% Plot the amplitude of the combined signal
cf = 0;

for i = 1 : length( SNR_simulated )
    signal_gen( SNR_simulated( i ) );
    [ SNR11, SNR12, SNR21, SNR22, rx_signal ] = decode();
    SNR_decoded = [ SNR_decoded ; [ SNR11, SNR12, SNR21, SNR22 ] ];
    
    cf = cf + 1;
    figure( cf ); 
    clf;
    plot( rx_signal );
    title( [ 'Rx Signal Amplitude (Simulated SNR: ' num2str( SNR_simulated( i ) ) ' dB)' ] );
end

% Plot the amplitude of the signal
cf = cf + 1;
figure( cf ); 
bh = bar( SNR_simulated , SNR_decoded );
grid on;
title( 'SNR decoded' );
xlabel( 'Simulated SNR (dB)' );
ylabel( 'Actual SNR (dB)' );