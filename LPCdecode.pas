{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit LPCdecode;

{$MODE Delphi}

interface
uses
iLBC_define,
helpfun,
lsf,
constants,C2Delphi_header;

   {---------------------------------------------------------------*
    *  interpolation of lsf coefficients for the decoder
    *--------------------------------------------------------------}
   procedure LSFinterpolate2a_dec(
       a:pareal;           { (o) lpc coefficients for a sub-frame }
       lsf1:pareal;    { (i) first lsf coefficient vector }
       lsf2:pareal;    { (i) second lsf coefficient vector }
       coef:real;         { (i) interpolation weight }
       length:integer          { (i) length of lsf vectors }
   );
   procedure SimplelsfDEQ(
       lsfdeq:pareal;    { (o) dequantized lsf coefficients }
       index:painteger;         { (i) quantization index }
       lpc_n:integer           { (i) number of LPCs }
   );
   procedure DecoderInterpolateLSF(
       syntdenum:pareal; { (o) synthesis filter coefficients }
       weightdenum:pareal; { (o) weighting denumerator
                                  coefficients }
       lsfdeq:pareal;       { (i) dequantized lsf coefficients }
       length:integer;         { (i) length of lsf coefficient vector }
       iLBCdec_inst:piLBC_Dec_Inst_t
                           { (i) the decoder state structure }
   );
implementation
   procedure LSFinterpolate2a_dec(
       a:pareal;           { (o) lpc coefficients for a sub-frame }
       lsf1:pareal;    { (i) first lsf coefficient vector }
       lsf2:pareal;    { (i) second lsf coefficient vector }
       coef:real;         { (i) interpolation weight }
       length:integer          { (i) length of lsf vectors }
   );
   var
   	 lsftmp:array [0..LPC_FILTERORDER] of real;
   begin
       interpolate(@lsftmp, lsf1, lsf2, coef, length);
       lsf2a(a, @lsftmp);
   end;

   {---------------------------------------------------------------*
    *  obtain dequantized lsf coefficients from quantization index
    *--------------------------------------------------------------}

   procedure SimplelsfDEQ(
       lsfdeq:pareal;    { (o) dequantized lsf coefficients }
       index:painteger;         { (i) quantization index }
       lpc_n:integer           { (i) number of LPCs }
   );
   var
       i, j, pos, cb_pos:integer;
   begin
       { decode first LSF }

       pos := 0;
       cb_pos := 0;
       for i := 0 to LSF_NSPLIT-1 do
       begin
           for j := 0 to dim_lsfCbTbl[i]-1 do
           begin
               lsfdeq[pos + j] := lsfCbTbl[cb_pos +(index[i])*dim_lsfCbTbl[i] + j];
           end;
           pos :=pos + dim_lsfCbTbl[i];
           cb_pos :=cb_pos + size_lsfCbTbl[i]*dim_lsfCbTbl[i];
       end;

       if (lpc_n>1) then
       begin
           { decode last LSF }

           pos := 0;
           cb_pos := 0;
           for i := 0 to LSF_NSPLIT-1 do
           begin
               for j := 0 to dim_lsfCbTbl[i]-1 do
               begin
                   lsfdeq[LPC_FILTERORDER + pos + j] := lsfCbTbl[cb_pos +(index[LSF_NSPLIT + i])*dim_lsfCbTbl[i] + j];
               end;
               pos :=pos + dim_lsfCbTbl[i];
               cb_pos :=cb_pos + size_lsfCbTbl[i]*dim_lsfCbTbl[i];
           end;
       end;
   end;

   {----------------------------------------------------------------*
    *  obtain synthesis and weighting filters form lsf coefficients
    *---------------------------------------------------------------}

   procedure DecoderInterpolateLSF(
       syntdenum:pareal; { (o) synthesis filter coefficients }
       weightdenum:pareal; { (o) weighting denumerator
                                  coefficients }
       lsfdeq:pareal;       { (i) dequantized lsf coefficients }
       length:integer;         { (i) length of lsf coefficient vector }
       iLBCdec_inst:piLBC_Dec_Inst_t
                           { (i) the decoder state structure }
   );
   var
       i, pos, lp_length:integer;
       lp:array [0..LPC_FILTERORDER ] of real;
       lsfdeq2:pareal;
   begin
       lsfdeq2 := @lsfdeq [ length];
       lp_length := length + 1;

       if (iLBCdec_inst^.mode=30) then
       begin
           { sub-frame 1: Interpolation between old and first }

           LSFinterpolate2a_dec(@lp, @iLBCdec_inst^.lsfdeqold, lsfdeq,
               lsf_weightTbl_30ms[0], length);
           move (lp,syntdenum[0],lp_length*sizeof(real));
           bwexpand(@weightdenum[0], @lp, LPC_CHIRP_WEIGHTDENUM,
               lp_length);

           { sub-frames 2 to 6: interpolation between first
              and last LSF }

           pos := lp_length;
           for i := 1 to 5 do
           begin
               LSFinterpolate2a_dec(@lp, lsfdeq, lsfdeq2,
                   lsf_weightTbl_30ms[i], length);
               move(lp,syntdenum[pos],lp_length*sizeof(real));
               bwexpand(@weightdenum [pos], @lp,
                   LPC_CHIRP_WEIGHTDENUM, lp_length);
               pos :=pos + lp_length;
           end;
       end
       else
       begin
           pos := 0;
           for i := 0 to iLBCdec_inst^.nsub-1 do
           begin
               LSFinterpolate2a_dec(@lp, @iLBCdec_inst^.lsfdeqold,
                   lsfdeq, lsf_weightTbl_20ms[i], length);
               move(lp,syntdenum[pos],lp_length*sizeof(real));
               bwexpand(@weightdenum[pos], @lp, LPC_CHIRP_WEIGHTDENUM,
                   lp_length);
               pos :=pos + lp_length;
           end;
       end;

       { update memory }

       if (iLBCdec_inst^.mode=30) then
           move( lsfdeq2[0],iLBCdec_inst^.lsfdeqold,length*sizeof(real))
       else
           move( lsfdeq[0],iLBCdec_inst^.lsfdeqold,length*sizeof(real));
   end;

end.
