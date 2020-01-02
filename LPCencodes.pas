unit LPCencodes;

{$MODE Delphi}

interface
uses
iLBC_define,
helpfun,
lsf,
constants,C2Delphi_header;

   {----------------------------------------------------------------*
    *  lpc analysis (subrutine to LPCencode)
    *---------------------------------------------------------------}
   procedure SimpleAnalysis(
       lsf:pareal;         { (o) lsf coefficients }
       data:pareal;    { (i) new data vector }
       iLBCenc_inst:piLBC_Enc_Inst_t
                           { (i/o) the encoder state structure }
   );
   procedure LSFinterpolate2a_enc(
       a:pareal;       { (o) lpc coefficients }
       lsf1:pareal;{ (i) first set of lsf coefficients }
       lsf2:pareal;{ (i) second set of lsf coefficients }
       coef:real;     { (i) weighting coefficient to use between
                              lsf1 and lsf2 }
       length:integer      { (i) length of coefficient vectors }
   );
   procedure SimpleInterpolateLSF(
       syntdenum:pareal;   { (o) the synthesis filter denominator
                                  resulting from the quantized
                                  interpolated lsf }
       weightdenum:pareal; { (o) the weighting filter denominator
                                  resulting from the unquantized
                                  interpolated lsf }
       lsf:pareal;         { (i) the unquantized lsf coefficients }
       lsfdeq:pareal;      { (i) the dequantized lsf coefficients }
       lsfold:pareal;      { (i) the unquantized lsf coefficients of
                                  the previous signal frame }
       lsfdeqold:pareal; { (i) the dequantized lsf coefficients of
                                  the previous signal frame }
       length:integer;         { (i) should equate LPC_FILTERORDER }
       iLBCenc_inst:piLBC_Enc_Inst_t
                           { (i/o) the encoder state structure }
   );
   procedure SimplelsfQ(
       lsfdeq:pareal;    { (o) dequantized lsf coefficients
                              (dimension FILTERORDER) }
       index:painteger;     { (o) quantization index }
       lsf:pareal;      { (i) the lsf coefficient vector to be
                              quantized (dimension FILTERORDER ) }
       lpc_n:integer     { (i) number of lsf sets to quantize }
   );
   procedure LPCencode(
       syntdenum:pareal; { (i/o) synthesis filter coefficients
                                  before/after encoding }
       weightdenum:pareal; { (i/o) weighting denumerator
                                  coefficients before/after
                                  encoding }
       lsf_index:painteger;     { (o) lsf quantization index }
       data:pareal;    { (i) lsf coefficients to quantize }
       iLBCenc_inst:piLBC_Enc_Inst_t
                           { (i/o) the encoder state structure }
   );
implementation
   procedure SimpleAnalysis(
       lsf:pareal;         { (o) lsf coefficients }
       data:pareal;    { (i) new data vector }
       iLBCenc_inst:piLBC_Enc_Inst_t
                           { (i/o) the encoder state structure }
   );
   var
       k, iss:integer;
       temp:array [0..BLOCKL_MAX-1] of real;
       lp:array [0..LPC_FILTERORDER] of real;
       lp2:array [0..LPC_FILTERORDER ] of real;
       r:array [0..LPC_FILTERORDER] of real;
   begin

       iss:=LPC_LOOKBACK+BLOCKL_MAX-iLBCenc_inst^.blockl;
       move(data[0],iLBCenc_inst^.lpc_buffer[iss],iLBCenc_inst^.blockl*sizeof(real));

       { No lookahead, last window is asymmetric }

       for k := 0 to iLBCenc_inst^.lpc_n-1 do
       begin
           iss := LPC_LOOKBACK;

           if (k < (iLBCenc_inst^.lpc_n - 1)) then
           begin
               window(@temp, @lpc_winTbl,@iLBCenc_inst^.lpc_buffer, BLOCKL_MAX);
           end
           else
           begin
               window(@temp, @lpc_asymwinTbl,@iLBCenc_inst^.lpc_buffer[iss], BLOCKL_MAX);
           end;
           autocorr(@r, @temp, BLOCKL_MAX, LPC_FILTERORDER);//這會爆
           window(@r, @r, @lpc_lagwinTbl, LPC_FILTERORDER + 1);
           levdurb(@lp, @temp, @r, LPC_FILTERORDER);
           bwexpand(@lp2, @lp, LPC_CHIRP_SYNTDENUM, LPC_FILTERORDER+1);
           a2lsf(@lsf [ k*LPC_FILTERORDER], @lp2);
       end;
       iss:=LPC_LOOKBACK+BLOCKL_MAX-iLBCenc_inst^.blockl;
       move(iLBCenc_inst^.lpc_buffer[LPC_LOOKBACK+BLOCKL_MAX-iss],iLBCenc_inst^.lpc_buffer[0],iss*sizeof(real));
   end;

   {----------------------------------------------------------------*
    *  lsf interpolator and conversion from lsf to a coefficients
    *  (subrutine to SimpleInterpolateLSF)
    *---------------------------------------------------------------}

   procedure LSFinterpolate2a_enc(
       a:pareal;       { (o) lpc coefficients }
       lsf1:pareal;{ (i) first set of lsf coefficients }
       lsf2:pareal;{ (i) second set of lsf coefficients }
       coef:real;     { (i) weighting coefficient to use between
                              lsf1 and lsf2 }
       length:integer      { (i) length of coefficient vectors }
   );
   var
   	lsftmp:array [0..LPC_FILTERORDER] of real;
   begin
       interpolate(@lsftmp, lsf1, lsf2, coef, length);
       lsf2a(a, @lsftmp);
   end;

   {----------------------------------------------------------------*
    *  lsf interpolator (subrutine to LPCencode)
    *---------------------------------------------------------------}

   procedure SimpleInterpolateLSF(
       syntdenum:pareal;   { (o) the synthesis filter denominator
                                  resulting from the quantized
                                  interpolated lsf }
       weightdenum:pareal; { (o) the weighting filter denominator
                                  resulting from the unquantized
                                  interpolated lsf }
       lsf:pareal;         { (i) the unquantized lsf coefficients }
       lsfdeq:pareal;      { (i) the dequantized lsf coefficients }
       lsfold:pareal;      { (i) the unquantized lsf coefficients of
                                  the previous signal frame }
       lsfdeqold:pareal; { (i) the dequantized lsf coefficients of
                                  the previous signal frame }
       length:integer;         { (i) should equate LPC_FILTERORDER }
       iLBCenc_inst:piLBC_Enc_Inst_t
                           { (i/o) the encoder state structure }
   );
   var
       i, pos, lp_length:integer;
       lp:array [0..LPC_FILTERORDER] of real;
       lsf2, lsfdeq2:pareal;
   begin

       lsf2 := @lsf [ length];
       lsfdeq2 := @lsfdeq [ length];
       lp_length := length + 1;

       if (iLBCenc_inst^.mode=30) then
       begin
           { sub-frame 1: Interpolation between old and first
              set of lsf coefficients }

           LSFinterpolate2a_enc(@lp[0], lsfdeqold, lsfdeq,
               lsf_weightTbl_30ms[0], length);
           move(lp[0],syntdenum[0],lp_length*sizeof(real));
           LSFinterpolate2a_enc(@lp[0], lsfold, lsf,
               lsf_weightTbl_30ms[0], length);
           bwexpand(@weightdenum[0], @lp[0], LPC_CHIRP_WEIGHTDENUM, lp_length);

           { sub-frame 2 to 6: Interpolation between first
              and second set of lsf coefficients }

           pos := lp_length;
           for i := 1 to iLBCenc_inst^.nsub-1 do
           begin
               LSFinterpolate2a_enc(@lp[0], lsfdeq, lsfdeq2,
                   lsf_weightTbl_30ms[i], length);
               move(lp[0],syntdenum [ pos],lp_length*sizeof(real));

               LSFinterpolate2a_enc(@lp[0], lsf, lsf2,
                   lsf_weightTbl_30ms[i], length);
               bwexpand(@weightdenum [ pos], @lp[0],
                   LPC_CHIRP_WEIGHTDENUM, lp_length);
               pos :=pos + lp_length;
           end;
       end
       else
       begin
           pos := 0;
           for i := 0 to iLBCenc_inst^.nsub-1 do
           begin
               LSFinterpolate2a_enc(@lp[0], lsfdeqold, lsfdeq,
                   lsf_weightTbl_20ms[i], length);
               move(lp[0],syntdenum[pos],lp_length*sizeof(real));
               LSFinterpolate2a_enc(@lp[0], lsfold, lsf,
                   lsf_weightTbl_20ms[i], length);
               bwexpand(@weightdenum[pos], @lp[0],
                   LPC_CHIRP_WEIGHTDENUM, lp_length);
               pos :=pos + lp_length;
           end;
       end;

       { update memory }

       if (iLBCenc_inst^.mode=30) then
       begin
           move( lsf2[0],lsfold[0], length*sizeof(real));
           move( lsfdeq2[0],lsfdeqold[0], length*sizeof(real));
       end
       else
       begin
           move( lsf[0],lsfold[0], length*sizeof(real));
           move( lsfdeq[0],lsfdeqold[0], length*sizeof(real));
       end;
   end;

   {----------------------------------------------------------------*
    *  lsf quantizer (subrutine to LPCencode)
    *---------------------------------------------------------------}

   procedure SimplelsfQ(
       lsfdeq:pareal;    { (o) dequantized lsf coefficients
                              (dimension FILTERORDER) }
       index:painteger;     { (o) quantization index }
       lsf:pareal;      { (i) the lsf coefficient vector to be
                              quantized (dimension FILTERORDER ) }
       lpc_n:integer     { (i) number of lsf sets to quantize }
   );
   begin
       { Quantize first LSF with memoryless split VQ }
       SplitVQ(lsfdeq, index, lsf, @lsfCbTbl, LSF_NSPLIT,
           @dim_lsfCbTbl, @size_lsfCbTbl);

       if (lpc_n=2) then
       begin
           { Quantize second LSF with memoryless split VQ }
           SplitVQ(@lsfdeq [ LPC_FILTERORDER], @index [ LSF_NSPLIT],
               @lsf [ LPC_FILTERORDER], @lsfCbTbl, LSF_NSPLIT,
               @dim_lsfCbTbl, @size_lsfCbTbl);
       end;
   end;

   {----------------------------------------------------------------*
    *  lpc encoder
    *---------------------------------------------------------------}

   procedure LPCencode(
       syntdenum:pareal; { (i/o) synthesis filter coefficients
                                  before/after encoding }
       weightdenum:pareal; { (i/o) weighting denumerator
                                  coefficients before/after
                                  encoding }
       lsf_index:painteger;     { (o) lsf quantization index }
       data:pareal;    { (i) lsf coefficients to quantize }
       iLBCenc_inst:piLBC_Enc_Inst_t
                           { (i/o) the encoder state structure }
   );
   var
       lsf:array [0..LPC_FILTERORDER * LPC_N_MAX-1] of real;
       lsfdeq:array [0..LPC_FILTERORDER * LPC_N_MAX-1] of real;
       change:integer;
   begin
   		//change:=0;
       SimpleAnalysis(@lsf, data, iLBCenc_inst);
       SimplelsfQ(@lsfdeq, lsf_index, @lsf, iLBCenc_inst^.lpc_n);
       change:=LSF_check(@lsfdeq, LPC_FILTERORDER, iLBCenc_inst^.lpc_n);
       SimpleInterpolateLSF(syntdenum, weightdenum,
           @lsf, @lsfdeq, @iLBCenc_inst^.lsfold,
           @iLBCenc_inst^.lsfdeqold, LPC_FILTERORDER, iLBCenc_inst);
   end;

end.