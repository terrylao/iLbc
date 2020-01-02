unit iLBC_define;

{$MODE Delphi}

interface

Const
   { general codec settings }

   FS             =         8000.0;
   BLOCKL_20MS    =         160          ;
   BLOCKL_30MS    =         240          ;
   BLOCKL_MAX     =         240          ;
   NSUB_20MS      =         4            ;
   NSUB_30MS      =         6            ;
   NSUB_MAX       =     6                ;
   NASUB_20MS     =         2            ;





   NASUB_30MS          =    4 ;
   NASUB_MAX           =    4 ;
   SUBL                =40    ;
   STATE_LEN           =    80;
   STATE_SHORT_LEN_30MS=    58;
   STATE_SHORT_LEN_20MS=    57;

   { LPC settings }

   LPC_FILTERORDER       =  10                 ;
   LPC_CHIRP_SYNTDENUM   =  0.9025             ;
   LPC_CHIRP_WEIGHTDENUM =  0.4222             ;
   LPC_LOOKBACK          =  60                 ;
   LPC_N_20MS            =  1                  ;
   LPC_N_30MS            =  2                  ;
   LPC_N_MAX             =  2                  ;
   LPC_ASYMDIFF          =  20                 ;
   LPC_BW                =  60.0               ;
   LPC_WN                =  1.0001             ;
   LSF_NSPLIT            =  3                  ;
   LSF_NUMBER_OF_STEPS   =  4                  ;
   LPC_HALFORDER         =  (LPC_FILTERORDER div 2);

   { cb settings }

   CB_NSTAGES        =      3         ;
   CB_EXPAND         =      2         ;
   CB_MEML           =      147       ;
   CB_FILTERLEN      =  2*4           ;
   CB_HALFFILTERLEN  =  4             ;
   CB_RESRANGE       =      34        ;
   CB_MAXGAIN        =      1.3       ;

   { enhancer }

   ENH_BLOCKL      =        80;  { block length }
   ENH_BLOCKL_HALF =        (ENH_BLOCKL div 2);
   ENH_HL          =        3;   { 2*ENH_HL+1 is number blocks
                                  in said second sequence }
   ENH_SLOP        =    2;   { max difference estimated and
                                   correct pitch period }
   ENH_PLOCSL      =        20;  { pitch-estimates and pitch-
                                   locations buffer length }
   ENH_OVERHANG    =    2;
   ENH_UPS0        =    4;   { upsampling rate }
   ENH_FL0         =        3;   { 2*FLO+1 is the length of
                                   each filter }
   ENH_VECTL       =        (ENH_BLOCKL+2*ENH_FL0);





   ENH_CORRDIM        =     (2*ENH_SLOP+1);
   ENH_NBLOCKS        =     (BLOCKL_MAX div ENH_BLOCKL);
   ENH_NBLOCKS_EXTRA  =     5;
   ENH_NBLOCKS_TOT    =     8;   { ENH_NBLOCKS +
                                           ENH_NBLOCKS_EXTRA }
   ENH_BUFL           = (ENH_NBLOCKS_TOT)*ENH_BLOCKL;
   ENH_ALPHA0         =     0.05;

   { Down sampling }

   FILTERORDER_DS    =      7;
   DELAY_DS          =  3    ;
   FACTOR_DS         =      2;

   { bit stream defs }

   NO_OF_BYTES_20MS =   38   ;
   NO_OF_BYTES_30MS =   50   ;
   NO_OF_WORDS_20MS =   19   ;
   NO_OF_WORDS_30MS =   25   ;
   STATE_BITS       =       3;
   BYTE_LEN         =   8    ;
   ULP_CLASSES      =       3;

   { help parameters }

   FLOAT_MAX      =        1.0e37                ;
   EPS            =        2.220446049250313e-016;
   PI             =        3.14159265358979323846;
   MIN_SAMPLE     =        -32768                ;
   MAX_SAMPLE     =        32767                 ;
   TWO_PI         =        6.283185307           ;
   PI2            =        0.159154943           ;
   
Type
   { type definition encoder instance }
   piLBC_ULP_Inst_t=^iLBC_ULP_Inst_t;
   iLBC_ULP_Inst_t=record 
       lsf_bits: array [0..5,0..ULP_CLASSES+1] of Integer;
       start_bits: array [0..ULP_CLASSES+1] of Integer;
       startfirst_bits:array [0..ULP_CLASSES+1] of Integer;
       scale_bits:array [0..ULP_CLASSES+1] of integer;
       state_bits:array [0..ULP_CLASSES+1] of integer;
       extra_cb_index:array [0..CB_NSTAGES-1,0..ULP_CLASSES+1] of integer;
       extra_cb_gain:array [0..CB_NSTAGES-1,0..ULP_CLASSES+1] of integer;
       cb_index:array [0..NSUB_MAX-1,0..CB_NSTAGES-1,0..ULP_CLASSES+1] of integer;
       cb_gain:array [0..NSUB_MAX-1,0..CB_NSTAGES-1,0..ULP_CLASSES+1] of integer;
   end;

   { type definition encoder instance }




   piLBC_Enc_Inst_t=^iLBC_Enc_Inst_t;
   iLBC_Enc_Inst_t=record

       { flag for frame size mode }
       mode:integer;

       { basic parameters for different frame sizes }
       blockl:integer;
       nsub:integer;
       nasub:integer;
       no_of_bytes, no_of_words:integer;
       lpc_n:integer;
       state_short_len:integer;
       ULP_inst:^iLBC_ULP_Inst_t;

       { analysis filter state }
       anaMem:array [0..LPC_FILTERORDER-1] of real;

       { old lsf parameters for interpolation }
       lsfold:array [0..LPC_FILTERORDER-1] of real;
       lsfdeqold:array [0..LPC_FILTERORDER-1] of real;

       { signal buffer for LP analysis }
       lpc_buffer:array [0..LPC_LOOKBACK + BLOCKL_MAX-1] of real;

       { state of input HP filter }
       hpimem:array [0..3] of real;

   end;

   { type definition decoder instance }
   piLBC_Dec_Inst_t=^iLBC_Dec_Inst_t;
   iLBC_Dec_Inst_t=record

       { flag for frame size mode }
       mode:integer;

       { basic parameters for different frame sizes }
       blockl:integer;
       nsub:integer;
       nasub:integer;
       no_of_bytes, no_of_words:integer;
       lpc_n:integer;
       state_short_len:integer;
       ULP_inst:^iLBC_ULP_Inst_t;

       { synthesis filter state }
       syntMem:array [0..LPC_FILTERORDER-1] of real;

       { old LSF for interpolation }


	   
       lsfdeqold:array [0..LPC_FILTERORDER-1 ] of real;

       { pitch lag estimated in enhancer and used in PLC }
       last_lag:integer;

       { PLC state information }
       prevLag, consPLICount, prevPLI, prev_enh_pl:integer;
       prevLpc:array [0..LPC_FILTERORDER] of real;
       prevResidual:array [0..NSUB_MAX*SUBL-1] of real;
       per:real;
       seed:cardinal;

       { previous synthesis filter parameters }
       old_syntdenum:array [0..(LPC_FILTERORDER + 1)*NSUB_MAX-1] of real;

       { state of output HP filter }
       hpomem:array [0..3] of real;

       { enhancer state information }
       use_enhancer:integer;
       enh_buf:array [0..ENH_BUFL-1] of real;
       enh_period:array [0..ENH_NBLOCKS_TOT-1] of real;

   end;

implementation
end.