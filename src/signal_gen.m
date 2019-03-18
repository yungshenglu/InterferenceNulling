function [ output ] = signal_gen( SNR_simulated )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wl_example_siso_ofdm_txrx.m
    % A detailed write-up of this example is available on the wiki:
    % http://warpproject.org/trac/wiki/WARPLab/Examples/OFDM
    %
    % Copyright (c) 2015 Mango Communications - All Rights Reserved
    % Distributed under the WARP License (http://warpproject.org/license)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Tx side
    % Params:
    % USE_WARPLAB_TXRX: Enable WARPLab-in-the-loop (otherwise sim-only)
    % CHANNEL: Channel to tune Tx and Rx radios
    USE_WARPLAB_TXRX        = 0;           
    CHANNEL                 = 11;          

    % Waveform params
    % N_OFDM_SYMS: Number of OFDM symbols
    % MOD_ORDER: Modulation order (2/4/16/64 = BSPK/QPSK/16-QAM/64-QAM)
    % TX_SCALE: Scale for Tx waveform ([ 0 : 1 ]) (TX gain)
    N_OFDM_SYMS             = 50;         
    MOD_ORDER               = 2;           
    TX_SCALE                = 1.0;         

    % OFDM params
    % SC_IND_PILOTS: Pilot subcarrier indices, index 1 is for DC and index 28 to 28 is not usesd form leakage.
    % SC_IND_DATA: Data subcarrier indices
    % N_SC: Number of subcarriers
    % CP_LEN: Cyclic prefix length
    % N_DATA_SYMS: Number of data symbols (one per data-bearing subcarrier per OFDM symbol)
    % INTERP_RATE: Interpolation rate (must be 2)
    SC_IND_PILOTS           = [ 8, 22, 44, 58 ];                                             
    SC_IND_DATA             = [ 2 : 7, 9 : 21, 23 : 27, 39 : 43, 45 : 57, 59 : 64 ];       
    N_SC                    = 64;                                                               
    CP_LEN                  = 16;                                                               
    N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);                                
    INTERP_RATE             = 2;                                                                

    % WARPLab experiment params
    % USE_AGC: Use the AGC if running on WARP hardware
    % MAX_TX_LEN: Maximum number of samples to use for this experiment
    % TRIGGER_OFFSET_TOL_NS: Trigger time offset toleration between Tx and Rx that can be accomodated
    USE_AGC                 = true;       
    MAX_TX_LEN              = 2 ^ 20;      
    TRIGGER_OFFSET_TOL_NS   = 3000;        

    if USE_WARPLAB_TXRX
        %% Set up the WARPLab experiment
        NUMNODES = 2;

        % Create a vector of node objects
        nodes   = wl_initNodes( NUMNODES );
        node_tx = nodes( 1 );
        node_rx = nodes( 2 );

        % Create a UDP broadcast trigger and tell each node to be ready for it
        eth_trig = wl_trigger_eth_udp_broadcast;
        wl_triggerManagerCmd( nodes, 'add_ethernet_trigger', [ eth_trig ] );

        % Read Trigger IDs into workspace
        trig_in_ids  = wl_getTriggerInputIDs( nodes( 1 ) );
        trig_out_ids = wl_getTriggerOutputIDs( nodes( 1 ) );

        % For both nodes, we will allow Ethernet to trigger the buffer baseband and the AGC
        wl_triggerManagerCmd( nodes, 'output_config_input_selection', [ trig_out_ids.BASEBAND, trig_out_ids.AGC ], [ trig_in_ids.ETH_A ] );

        % Set the trigger output delays.
        nodes.wl_triggerManagerCmd( 'output_config_delay', [ trig_out_ids.BASEBAND], 0 );
        nodes.wl_triggerManagerCmd( 'output_config_delay', [ trig_out_ids.AGC ], TRIGGER_OFFSET_TOL_NS );

        % Get IDs for the interfaces on the boards. 
        ifc_ids_TX = wl_getInterfaceIDs( node_tx );
        ifc_ids_RX = wl_getInterfaceIDs( node_rx );

        % Set up the TX / RX nodes and RF interfaces
        TX_RF     = ifc_ids_TX.RF_A;
        TX_RF_VEC = ifc_ids_TX.RF_A;
        TX_RF_ALL = ifc_ids_TX.RF_ALL;

        RX_RF     = ifc_ids_RX.RF_A;
        RX_RF_VEC = ifc_ids_RX.RF_A;
        RX_RF_ALL = ifc_ids_RX.RF_ALL;

        % Set up the interface for the experiment
        wl_interfaceCmd( node_tx, TX_RF_ALL , 'channel', 2.4, CHANNEL );
        wl_interfaceCmd( node_rx, RX_RF_ALL , 'channel', 2.4, CHANNEL );
        wl_interfaceCmd( node_tx, TX_RF_ALL , 'tx_gains', 3, 30 );

        if USE_AGC
            wl_interfaceCmd( node_rx, RX_RF_ALL, 'rx_gain_mode', 'automatic' );
            wl_basebandCmd( nodes, 'agc_target', -13 );
        else
            wl_interfaceCmd( node_rx, RX_RF_ALL, 'rx_gain_mode', 'manual' );
            % Rx RF Gain in [ 1 : 3 ]
            RxGainRF = 2;    
            % Rx Baseband Gain in [ 0 : 31 ]
            RxGainBB = 12;                 
            wl_interfaceCmd( node_rx, RX_RF_ALL, 'rx_gains', RxGainRF, RxGainBB );
        end

        % Get parameters from the node
        SAMP_FREQ    = wl_basebandCmd( nodes( 1 ), 'tx_buff_clk_freq' );
        Ts           = 1 / SAMP_FREQ;

        % We will read the transmitter's maximum I/Q buffer length and assign that value to a temporary variable.
        % NOTE:  We assume that the buffers sizes are the same for all interfaces

        maximum_buffer_len = min( MAX_TX_LEN, wl_basebandCmd( node_tx, TX_RF_VEC, 'tx_buff_max_num_samples' ) );
        example_mode_string = 'hw';
    else
        % Use same defaults for hardware-dependent params in sim-only version
        maximum_buffer_len  = min( MAX_TX_LEN, 2 ^ 20 );
        SAMP_FREQ           = 40e6;
        example_mode_string = 'sim';
    end

    
    %% David: Random two channel values for h1 and h2
    h1 = ( randn( 1, N_SC ) + 1i * randn( 1, N_SC ) );
    h1 = h1 ./ abs( h1 );
    h2 = ( randn( 1, N_SC ) + 1i * randn( 1, N_SC ) );
    h2 = h2 ./ abs( h2 );
    % end
    
    
    %% David: Calculate two precoding coefficients for Tx1 and Tx2
    w1 = ones( 1, N_SC );
    w2 = -( h1./ h2 ) .* w1;
    w_power = abs( w1 ) .^ 2 + abs( w2 ) .^ 2;
    w1 = w1 ./ sqrt( w_power );
    w2 = w2 ./ sqrt( w_power );
    % end

    
    %% Define a half-band 2x interpolation filter response
    interp_filt2 = zeros( 1, 43 );
    interp_filt2( [ 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21 ] ) = [ 12, -32, 72, -140, 252, -422, 682, -1086, 1778, -3284, 10364 ];
    interp_filt2( [ 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43 ] ) = interp_filt2( fliplr( [ 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21 ] ) );
    interp_filt2( 22 ) = 16384;
    interp_filt2 = interp_filt2 ./ max( abs( interp_filt2 ) );

    
    %% Define the preamble
    % Note: The STS symbols (short preamble) in the preamble meet the requirements needed by the
    % AGC core at the receiver. Details on the operation of the AGC are
    % available on the wiki: http://warpproject.org/trac/wiki/WARPLab/AGC
    sts_f = zeros( 1, 64 );
    sts_f( 1 : 27 ) = [ 0, 0, 0, 0, -1-1i, 0, 0, 0, -1-1i, 0, 0, 0, 1+1i, 0, 0, 0, 1+1i, 0, 0, 0, 1+1i, 0, 0, 0, 1+1i, 0, 0 ];
    sts_f( 39 : 64 ) = [ 0, 0, 1+1i, 0, 0, 0, -1-1i, 0, 0, 0, 1+1i, 0, 0, 0, -1-1i, 0, 0, 0, -1-1i, 0, 0, 0, 1+1i, 0, 0, 0 ];
    
    
    %% David: Multiply short preamble by each channel value
    sts1_t = ifft( sqrt( 13 / 6 ) .* sts_f .* h1, 64 );
    sts1_t = sts1_t( 1 : 16 );
    sts2_t = ifft( sqrt( 13 / 6 ) .* sts_f .* h2, 64 );
    sts2_t = sts1_t( 1 : 16 );
    % end

    %% LTS for CFO and channel estimation
    lts_f = [ 0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1 ];
    
    
    %% David: Multiply long preamble by each channel value
    lts1_t = lts_f .* h1;
    lts2_t = lts_f .* h2;
    lts1_t = ifft( lts1_t, 64 );
    lts2_t = ifft( lts2_t, 64 );
    % end

    
    %% David: Use 30 copies of the 16-sample STS for extra AGC settling margin for each Tx
    preamble1 = [ repmat( sts1_t, 1, 30 ),  lts1_t( 33 : 64 ), lts1_t, lts1_t ];
    preamble2 = [ repmat( sts2_t, 1, 30 ),  lts2_t( 33 : 64 ), lts2_t, lts2_t ];
    % end
    

    %% Sanity check variables that affect the number of Tx samples
    num_samps_needed = ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) + ...
                        INTERP_RATE * ( ( N_OFDM_SYMS * ( N_SC + CP_LEN ) ) + length( preamble1 ) +  ceil( length( interp_filt2 ) / 2 ) );

    if num_samps_needed > maximum_buffer_len
        fprintf( 'Too many OFDM symbols for TX_NUM_SAMPS!\n' );
        fprintf( 'Raise MAX_TX_LEN to %d, or \n', num_samps_needed );
        fprintf( 'Reduce N_OFDM_SYMS to %d\n', floor( ( ( maximum_buffer_len - ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) ) / INTERP_RATE - ( length( preamble ) +  ceil( length( interp_filt2 ) / 2 ) ) ) / ( N_SC + CP_LEN ) ) );
        return;
    end


    %% Generate a payload of random integers
    tx_data = randi( MOD_ORDER, 1, N_DATA_SYMS ) - 1;

    % Functions for data -> complex symbol mapping (like qammod, avoids comm toolbox requirement)
    % These anonymous functions implement the modulation mapping from IEEE 802.11-2012 Section 18.3.5.8
    modvec_bpsk   =  1 .* [ -1, 1 ];
    modvec_16qam  =  ( 1 / sqrt( 10 ) ) .* [ -3, -1, +3, +1 ];
    modvec_64qam  =  ( 1 / sqrt( 43 ) ) .* [ -7, -5, -1, -3, +7, +5, +1, +3 ];

    mod_fcn_bpsk  = @( x ) complex( modvec_bpsk( 1 + x ), 0 );
    mod_fcn_qpsk  = @( x ) complex( modvec_bpsk( 1 + bitshift( x, -1 ) ), modvec_bpsk( 1 + mod( x, 2 ) ) );
    mod_fcn_16qam = @( x ) complex( modvec_16qam( 1 + bitshift( x, -2 ) ), modvec_16qam( 1 + mod( x, 4 ) ) );
    mod_fcn_64qam = @( x ) complex( modvec_64qam( 1 + bitshift( x, -3 ) ), modvec_64qam( 1 + mod( x, 8 ) ) );

    %% Map the data values on to complex symbols
    switch MOD_ORDER
        % BPSK
        case 2         
            tx_syms = arrayfun( mod_fcn_bpsk, tx_data );
        % QPSK
        case 4         
            tx_syms = arrayfun( mod_fcn_qpsk, tx_data );
        % 16-QAM
        case 16        
            tx_syms = arrayfun( mod_fcn_16qam, tx_data );
        % 64-QAM
        case 64        
            tx_syms = arrayfun( mod_fcn_64qam, tx_data );
        otherwise
            fprintf( 'Invalid MOD_ORDER (%d)!  Must be in [2, 4, 16, 64]\n', MOD_ORDER );
            return;
    end

    % Reshape the symbol vector to a matrix with one column per OFDM symbol
    tx_syms_mat = reshape( tx_syms, length( SC_IND_DATA ), N_OFDM_SYMS );

    % Define the pilot tone values as BPSK symbols
    pilots = [ 1, 1, -1, 1 ].';

    % Repeat the pilots across all OFDM symbols
    pilots_mat = repmat( pilots, 1, N_OFDM_SYMS );


    %% IFFT
    % Construct the IFFT input matrix
    ifft_in_mat = zeros( N_SC, N_OFDM_SYMS );

    % Insert the data and pilot values; other subcarriers will remain at 0
    ifft_in_mat( SC_IND_DATA, : ) = tx_syms_mat;
    ifft_in_mat( SC_IND_PILOTS, : ) = pilots_mat;

    
    %% David: Perform the IFFT with precoding and multiply by each channel value
    tx_payload_w_mat1 = ifft( ifft_in_mat .* h1.' .* w1.', N_SC, 1 );
    tx_payload_w_mat2 = ifft( ifft_in_mat .* h2.' .* w2.', N_SC, 1 );
    % end
    
    
    %% David: Perform the IFFT without precoding and multiply by each channel value
    tx_payload_wo_mat1 = ifft( ifft_in_mat .* h1.', N_SC, 1 );
    tx_payload_wo_mat2 = ifft( ifft_in_mat .* h2.', N_SC, 1 );
    % end

    
    %% David: Insert the cyclic prefix
    if CP_LEN > 0
        tx_w_cp1 = tx_payload_w_mat1( ( end - CP_LEN + 1 : end ), : );
        tx_payload_w_mat1 = [ tx_w_cp1 ; tx_payload_w_mat1 ];
        
        tx_w_cp2 = tx_payload_w_mat2( ( end - CP_LEN + 1 : end ), : );
        tx_payload_w_mat2 = [ tx_w_cp2 ; tx_payload_w_mat2 ];
        
        tx_wo_cp1 = tx_payload_wo_mat1( ( end - CP_LEN + 1 : end ), : );
        tx_payload_wo_mat1 = [ tx_wo_cp1 ; tx_payload_wo_mat1 ];
        
        tx_wo_cp2 = tx_payload_wo_mat2( ( end - CP_LEN + 1 : end ), : );
        tx_payload_wo_mat2 = [ tx_wo_cp2 ; tx_payload_wo_mat2 ];
    end
    % end
    

    %% David: Reshape to a vector
    tx_payload_w_vec1 = reshape( tx_payload_w_mat1, 1, numel( tx_payload_w_mat1 ) );
    tx_payload_w_vec2 = reshape( tx_payload_w_mat2, 1, numel( tx_payload_w_mat2 ) );
    tx_payload_wo_vec1 = reshape( tx_payload_wo_mat1, 1, numel( tx_payload_wo_mat1 ) );
    tx_payload_wo_vec2 = reshape( tx_payload_wo_mat2, 1, numel( tx_payload_wo_mat2 ) );
    % end

    
    %% David: Construct the full time-domain OFDM waveform
    tx_vec1 = [ preamble1, zeros( 1, length( preamble2 ) ), tx_payload_wo_vec1, zeros( 1, length( tx_payload_wo_vec2 ) ), tx_payload_wo_vec1 / sqrt( 2 ), tx_payload_w_vec1 ];
    tx_vec2 = [ zeros( 1, length( preamble1 ) ), preamble2, zeros( 1, length( tx_payload_wo_vec1 ) ), tx_payload_wo_vec2, tx_payload_wo_vec2 / sqrt( 2 ), tx_payload_w_vec2 ];
    % end
    

    %% David: No need in this lab
    %{
    % Pad with zeros for transmission to deal with delay through the interpolation filter
    tx_vec_padded = [ tx_vec, zeros( 1, ceil( length( interp_filt2 ) / 2 ) ) ];


    %% Interpolate
    % Zero pad then filter (same as interp or upfirdn without signal processing toolbox)
    if INTERP_RATE ~= 2
       fprintf( 'Error: INTERP_RATE must equal 2\n' ); 
       return;
    end

    tx_vec_2x = zeros( 1 , 2 * numel( tx_vec_padded ) );
    tx_vec_2x( 1 : 2 : end ) = tx_vec_padded;
    tx_vec_air = filter( interp_filt2 , 1 , tx_vec_2x );
    %}
    % end
    
    
    %% David: Scale the Tx1 and Tx2 vector to +/- 1
    tx_vec_air1 = TX_SCALE .* tx_vec1 ./ max( abs( [ tx_vec1, tx_vec2 ] ) );
    tx_vec_air2 = TX_SCALE .* tx_vec2 ./ max( abs( [ tx_vec1, tx_vec2 ] ) );
    tx_vec_air1 = [ tx_vec_air1, zeros( 1, ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) ) ];
    tx_vec_air2 = [ tx_vec_air2, zeros( 1, ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) ) ];
    % end

    
    TX_NUM_SAMPS = length( tx_vec_air1 );
    if USE_WARPLAB_TXRX
        wl_basebandCmd( nodes, 'tx_delay', 0 );
        % Number of samples to send
        wl_basebandCmd( nodes, 'tx_length', TX_NUM_SAMPS );
        % Number of samples to receive
        wl_basebandCmd( nodes, 'rx_length', TX_NUM_SAMPS + ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) );
    end

    if USE_WARPLAB_TXRX
        % Write the Tx waveform to the Tx node
        wl_basebandCmd( node_tx, TX_RF_VEC, 'write_IQ', tx_vec_air( : ) );

        % Enable the Tx and Rx radios
        wl_interfaceCmd( node_tx, TX_RF, 'tx_en' );
        wl_interfaceCmd( node_rx, RX_RF, 'rx_en' );

        % Enable the Tx and Rx buffers
        wl_basebandCmd( node_tx, TX_RF, 'tx_buff_en' );
        wl_basebandCmd( node_rx, RX_RF, 'rx_buff_en' );

        % Trigger the Tx/Rx cycle at both nodes
        eth_trig.send();

        % Retrieve the received waveform from the Rx node
        rx_vec_air = wl_basebandCmd( node_rx, RX_RF_VEC, 'read_IQ', 0, TX_NUM_SAMPS + ( ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) ) );

        rx_vec_air = rx_vec_air( : ).';

        % Disable the Tx/Rx radios and buffers
        wl_basebandCmd( node_tx, TX_RF_ALL, 'tx_rx_buff_dis' );
        wl_basebandCmd( node_rx, RX_RF_ALL, 'tx_rx_buff_dis' );

        wl_interfaceCmd( node_tx, TX_RF_ALL, 'tx_rx_dis' );
        wl_interfaceCmd( node_rx, RX_RF_ALL, 'tx_rx_dis' );
    else
        % Sim-only mode: Apply wireless degradations here for sim (noise, fading, etc)

        % Perfect (i.e., Rx = Tx1 + Tx2):
        rx_vec_air = tx_vec_air1 + tx_vec_air2;

        % AWGN:
        % rx_vec_air = [ tx_vec_air, zeros( 1, ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) ) ];

        
        %% David: Simulated SNR: 5, 10, 15, 20 and 25 dB.
        SNR_dB = SNR_simulated;
        SNR = 10 ^ ( SNR_dB / 10 );
        Pn = ( mean( abs( tx_vec_air1 ) .^ 2 + abs( tx_vec_air2 ) .^ 2 ) ) / SNR;
        rx_vec_air = rx_vec_air + sqrt( Pn ) * complex( randn( 1, length( rx_vec_air ) ), randn( 1, length( rx_vec_air ) ) ) / sqrt( 2 );
        % end
        
        
        % CFO:
        % rx_vec_air = tx_vec_air .* exp( -1i * 2 * pi * 1e-4 * [ 0 : length( tx_vec_air ) - 1 ] );
    end


    %% Output files
    file = fopen( '../out/tx_data.bin', 'w' );
    fwrite( file, tx_data, 'int' );
    fclose( file );

    z_real = real( tx_syms_mat );
    z_imag = imag( tx_syms_mat );
    adjacent = [ z_real, z_imag ];
    file = fopen( '../out/tx_syms_mat.bin', 'w' );
    fwrite( file, adjacent, 'float' );
    fclose( file );

    complex_array = zeros( 1, 2 * length( rx_vec_air ) );
    for j = 1 : length( rx_vec_air )
        complex_array( 2 * j - 1 ) = real( rx_vec_air( j ) );
        complex_array( 2 * j ) = imag( rx_vec_air( j ) );
    end
    file = fopen( '../out/rx_vec_air.bin', 'w' );
    fwrite( file, complex_array, 'float' );
    fclose( file );
end
