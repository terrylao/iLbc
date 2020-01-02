{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit StateSearchWs;

{$MODE Delphi}

interface
uses
iLBC_define,
constants,
filter,
helpfun,C2Delphi_header,Math;

   {----------------------------------------------------------------*
    *  predictive noise shaping encoding of scaled start state
    *  (subrutine for StateSearchW)
    *---------------------------------------------------------------}
   procedure AbsQuantW(
       iLBCenc_inst:piLBC_Enc_Inst_t;
                           { (i) Encoder instance }
       nin:pareal;          { (i) vector to encode }
       syntDenum:pareal;   { (i) denominator of synthesis filter }
       weightDenum:pareal; { (i) denominator of weighting filter }
       nout:painteger;           { (o) vector of quantizer indexes }
       len:integer;        { (i) length of vector to encode and
                                  vector of quantizer indexes }
       state_first:integer     { (i) position of start state in the
                                  80 vec }
   );
   procedure StateSearchW(
       iLBCenc_inst:piLBC_Enc_Inst_t;
                           { (i) Encoder instance }
       residual:pareal;{ (i) target residual vector }
       syntDenum:pareal;   { (i) lpc synthesis filter }
       weightDenum:pareal; { (i) weighting filter denuminator }
       idxForMax:pinteger;     { (o) quantizer index for maximum
                                  amplitude }
       idxVec:painteger;    { (o) vector of quantization indexes }
       len:integer;        { (i) length of all vectors }
       state_first:integer     { (i) position of start state in the
                                  80 vec }
   );
implementation
   procedure AbsQuantW(
       iLBCenc_inst:piLBC_Enc_Inst_t;
                           { (i) Encoder instance }
       nin:pareal;          { (i) vector to encode }
       syntDenum:pareal;   { (i) denominator of synthesis filter }
       weightDenum:pareal; { (i) denominator of weighting filter }
       nout:painteger;           { (o) vector of quantizer indexes }
       len:integer;        { (i) length of vector to encode and
                                  vector of quantizer indexes }
       state_first:integer     { (i) position of start state in the
                                  80 vec }
   );
   var
       syntOut:pareal;
       syntOutBuf:array [0..LPC_FILTERORDER+STATE_SHORT_LEN_30MS-1] of real;
       toQ, xq:real;
       n:integer;
       index:integer;
   begin

       { initialization of buffer for filtering }

       fillchar(syntOutBuf, LPC_FILTERORDER*sizeof(real), 0);

       { initialization of pointer for filtering }

       syntOut := @syntOutBuf[LPC_FILTERORDER];

       { synthesis and weighting filters on input }

       if (state_first<>0) then
       begin
           AllPoleFilter (nin, weightDenum, SUBL, LPC_FILTERORDER);
       end
       else
       begin
           AllPoleFilter (nin, weightDenum,iLBCenc_inst^.state_short_len-SUBL,LPC_FILTERORDER);
       end;

       { encoding loop }

       for n:=0 to len-1 do
       begin
           { time update of filter coefficients }

           if ((state_first<>0) and (n=SUBL))then
           begin
               syntDenum := @syntDenum [ (LPC_FILTERORDER+1)];
               weightDenum :=@ weightDenum [ (LPC_FILTERORDER+1)];

               { synthesis and weighting filters on input }
               AllPoleFilter (@nin[n], weightDenum, len-n,LPC_FILTERORDER);

           end
           else 
           if ((state_first=0) and (n=(iLBCenc_inst^.state_short_len-SUBL))) then
           begin
               syntDenum :=@syntDenum[ (LPC_FILTERORDER+1)];
               weightDenum :=@weightDenum[ (LPC_FILTERORDER+1)];

               { synthesis and weighting filters on input }
               AllPoleFilter (@nin[n], weightDenum, len-n, LPC_FILTERORDER);
           end;

           { prediction of synthesized and weighted input }

           syntOut[n] := 0.0;
           AllPoleFilter (@syntOut[n], weightDenum, 1,LPC_FILTERORDER);

           { quantization }

           toQ := nin[n]-syntOut[n];

           sort_sq(@xq, @index, toQ, @state_sq3Tbl, 8);
           nout[n]:=index;
           syntOut[n] := state_sq3Tbl[nout[n]];

           { update of the prediction filter }

           AllPoleFilter(@syntOut[n], weightDenum, 1,LPC_FILTERORDER);
       end;
   end;

   {----------------------------------------------------------------*
    *  encoding of start state
    *---------------------------------------------------------------}

   procedure StateSearchW(
       iLBCenc_inst:piLBC_Enc_Inst_t;
                           { (i) Encoder instance }
       residual:pareal;{ (i) target residual vector }
       syntDenum:pareal;   { (i) lpc synthesis filter }
       weightDenum:pareal; { (i) weighting filter denuminator }
       idxForMax:pinteger;     { (o) quantizer index for maximum
                                  amplitude }
       idxVec:painteger;    { (o) vector of quantization indexes }
       len:integer;        { (i) length of all vectors }
       state_first:integer     { (i) position of start state in the
                                  80 vec }
   );
   var
       dtmp, maxVal:real;
       tmpbuf:array [0..LPC_FILTERORDER+2*STATE_SHORT_LEN_30MS-1] of real;
       tmp:pareal;
       numerator:array [0..LPC_FILTERORDER] of real;
       foutbuf:array [0..LPC_FILTERORDER+2*STATE_SHORT_LEN_30MS-1] of real;
       fout:pareal;
       k:integer;
       qmax, scal:real;
   begin

       { initialization of buffers and filter coefficients }

       fillchar(tmpbuf, LPC_FILTERORDER*sizeof(real), 0);
       fillchar(foutbuf, LPC_FILTERORDER*sizeof(real), 0);
       for k:=0 to LPC_FILTERORDER-1 do
       begin
           numerator[k]:=syntDenum[LPC_FILTERORDER-k];
       end;
       numerator[LPC_FILTERORDER]:=syntDenum[0];
       tmp := @tmpbuf[LPC_FILTERORDER];
       fout := @foutbuf[LPC_FILTERORDER];

       { circular convolution with the all-pass filter }

       move( residual[0],tmp[0], len*sizeof(real));
       fillchar(tmp[len], len*sizeof(real), 0);
       ZeroPoleFilter(tmp, @numerator, syntDenum, 2*len,LPC_FILTERORDER, fout);
       for k:=0 to len-1 do
       begin
           fout[k] :=fout[k] + fout[k+len];
       end;

       { identification of the maximum amplitude value }

       maxVal := fout[0];
       for k:=1 to len-1 do
       begin
           if (fout[k]*fout[k] > maxVal*maxVal) then
           begin
               maxVal := fout[k];
           end;
       end;
       maxVal:=abs(maxVal);

       { encoding of the maximum amplitude value }

       if (maxVal < 10.0) then
       begin
           maxVal := 10.0;
       end;
       maxVal := log10(maxVal);
       sort_sq(@dtmp, idxForMax, maxVal, @state_frgqTbl, 64);

       { decoding of the maximum amplitude representation value,
          and corresponding scaling of start state }

       maxVal:=state_frgqTbl[idxForMax^];
       qmax := power(10,maxVal);
       scal := (4.5)/qmax;
       for k:=0 to len-1 do
       begin
           fout[k] :=fout[k] * scal;
       end;

       { predictive noise shaping encoding of scaled start state }

       AbsQuantW(iLBCenc_inst, fout,syntDenum, weightDenum,idxVec, len, state_first);
   end;
end.
