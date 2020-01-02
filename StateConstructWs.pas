unit StateConstructWs;

{$MODE Delphi}

interface
uses
iLBC_define,
constants,
filter,C2Delphi_header,math;

   {----------------------------------------------------------------*
    *  decoding of the start state
    *---------------------------------------------------------------}
   procedure StateConstructW(
       idxForMax:integer;      { (i) 6-bit index for the quantization of
                                  max amplitude }
       idxVec:painteger;    { (i) vector of quantization indexes }
       syntDenum:pareal;   { (i) synthesis filter denumerator }
       nout:pareal;         { (o) the decoded state vector }
       len:integer             { (i) length of a state vector }
   );
implementation
   procedure StateConstructW(
       idxForMax:integer;      { (i) 6-bit index for the quantization of
                                  max amplitude }
       idxVec:painteger;    { (i) vector of quantization indexes }
       syntDenum:pareal;   { (i) synthesis filter denumerator }
       nout:pareal;         { (o) the decoded state vector }
       len:integer             { (i) length of a state vector }
   );
   var
       maxVal:real;
       tmpbuf:array [0..LPC_FILTERORDER+2*STATE_LEN-1] of real;
       tmp:pareal;
       numerator:array [0..LPC_FILTERORDER] of real;
       foutbuf:array [0..LPC_FILTERORDER+2*STATE_LEN-1] of real;
       fout:pareal;
       k,tmpi:integer;
   begin

       { decoding of the maximum value }

       maxVal := state_frgqTbl[idxForMax];
       maxVal := power(10,maxVal)/4.5;

       { initialization of buffers and coefficients }

       fillchar(tmpbuf, LPC_FILTERORDER*sizeof(real), 0);
       fillchar(foutbuf, LPC_FILTERORDER*sizeof(real), 0);
       for k:=0 to LPC_FILTERORDER-1 do
       begin
           numerator[k]:=syntDenum[LPC_FILTERORDER-k];
       end;
       numerator[LPC_FILTERORDER]:=syntDenum[0];
       tmp := @tmpbuf[LPC_FILTERORDER];
       fout := @foutbuf[LPC_FILTERORDER];

       { decoding of the sample values }

       for k:=0 to len-1 do
       begin
           tmpi := len-1-k;
           { maxVal := 1/scal }
           tmp[k] := maxVal*state_sq3Tbl[idxVec[tmpi]];
       end;

       { circular convolution with all-pass filter }

       fillchar(tmp[len], len*sizeof(real), 0);
       ZeroPoleFilter(tmp, @numerator, syntDenum, 2*len,LPC_FILTERORDER, fout);
       for k:=0 to len-1 do
       begin
           nout[k] := fout[len-1-k]+fout[2*len-1-k];
       end;
   end;
end.