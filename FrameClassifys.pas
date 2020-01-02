{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit FrameClassifys;

{$MODE Delphi}

interface
uses iLBC_define,C2Delphi_header;

   {---------------------------------------------------------------*
    *  Classification of subframes to localize start state
    *--------------------------------------------------------------}
   function FrameClassify(      { index to the max-energy sub-frame }
       iLBCenc_inst:piLBC_Enc_Inst_t;
                           { (i/o) the encoder state structure }
       residual:pareal     { (i) lpc residual signal }
   ):integer;
implementation
   function FrameClassify(      { index to the max-energy sub-frame }
       iLBCenc_inst:piLBC_Enc_Inst_t;
                           { (i/o) the encoder state structure }
       residual:pareal     { (i) lpc residual signal }
   ):integer;
   var
       max_ssqEn:real;
       fssqEn:array [0..NSUB_MAX-1] of real;
       bssqEn:array [0..NSUB_MAX-1] of real;
       pp:^real;
       n, l, max_ssqEn_n:integer;
       ssqEn_win:array [0..NSUB_MAX-1-1] of real;
       sampEn_win:array [0..4] of real;
   begin
   		ssqEn_win[0]:=0.8;
   		ssqEn_win[1]:=0.9;
   		ssqEn_win[2]:=1.0;
   		ssqEn_win[3]:=0.9;
   		ssqEn_win[4]:=0.8;
      sampEn_win[0]:=1.0/6.0;
      sampEn_win[1]:=2.0/6.0;
      sampEn_win[2]:=3.0/6.0;
      sampEn_win[3]:=4.0/6.0;
      sampEn_win[4]:=5.0/6.0;

       { init the front and back energies to zero }

       fillchar(fssqEn, NSUB_MAX*sizeof(real), 0);
       fillchar(bssqEn, NSUB_MAX*sizeof(real), 0);

       { Calculate front of first seqence }

       n:=0;
       pp:=@residual[0];
       for l:=0 to 4 do
       begin
           fssqEn[n] :=fssqEn[n] + sampEn_win[l] * pp^ * pp^;
           inc(pp);
       end;
       for l:=5 to SUBL-1 do
       begin
           fssqEn[n] :=fssqEn[n]+ (pp^) * (pp^);
           inc(pp);
       end;

       { Calculate front and back of all middle sequences }

       for n:=1 to iLBCenc_inst^.nsub-2 do
       begin
           pp:=@residual[n*SUBL];
           for l:=0 to 4 do
           begin
               fssqEn[n] :=fssqEn[n] + sampEn_win[l] * (pp^) * (pp^);
               bssqEn[n] :=bssqEn[n] + (pp^) * (pp^);
               inc(pp);
           end;
           for l:=5 to SUBL-6 do
           begin
               fssqEn[n] :=fssqEn[n] + (pp^) * (pp^);
               bssqEn[n] :=bssqEn[n] + (pp^) * (pp^);
               inc(pp);
           end;
           for l:=SUBL-5 to SUBL-1 do
           begin
               fssqEn[n] :=fssqEn[n]+ (pp^) * (pp^);
               bssqEn[n] :=bssqEn[n]+ sampEn_win[SUBL-l-1] * (pp^) * (pp^);
               inc(pp);
           end;
       end;

       { Calculate back of last seqence }

       n:=iLBCenc_inst^.nsub-1;
       pp:=@residual[n*SUBL];
       for l:=0 to SUBL-6 do
       begin
           bssqEn[n] :=bssqEn[n] + (pp^) * (pp^);
           inc(pp);
       end;
       for l:=SUBL-5 to SUBL-1 do
       begin
           bssqEn[n] :=bssqEn[n] + sampEn_win[SUBL-l-1] * (pp^) * (pp^);
           inc(pp);
       end;

       { find the index to the weighted 80 sample with
          most energy }

       if (iLBCenc_inst^.mode=20) then
       	l:=1
       else
        l:=0;

       max_ssqEn:=(fssqEn[0]+bssqEn[1])*ssqEn_win[l];
       max_ssqEn_n:=1;
       for n:=2 to iLBCenc_inst^.nsub-1 do
       begin
           inc(l);
           if ((fssqEn[n-1]+bssqEn[n])*ssqEn_win[l] > max_ssqEn) then
           begin
               max_ssqEn:=(fssqEn[n-1]+bssqEn[n]) * ssqEn_win[l];
               max_ssqEn_n:=n;
           end;
       end;
       result := max_ssqEn_n;
   end;
end.
