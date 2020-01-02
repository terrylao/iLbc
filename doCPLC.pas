{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit doCPLC;

{$MODE Delphi}

interface
uses iLBC_define,C2Delphi_header;
   {----------------------------------------------------------------*
    *  Compute cross correlation and pitch gain for pitch prediction
    *  of last subframe at given lag.
    *---------------------------------------------------------------}
   procedure compCorr(
       cc:pareal;      { (o) cross correlation coefficient }
       gc:pareal;      { (o) gain }
       pm:pareal;
       buffer:PAreal;  { (i) signal buffer }
       lag:integer;    { (i) pitch lag }
       bLen:integer;       { (i) length of buffer }
       sRange:integer      { (i) correlation search length }
   );
   procedure doThePLC(
       PLCresidual:pareal; { (o) concealed residual }
       PLClpc:pareal;      { (o) concealed LP parameters }
       PLI:integer;        { (i) packet loss indicator
                                  0 - no PL, 1 = PL }
       decresidual:pareal; { (i) decoded residual }
       lpc:pareal;         { (i) decoded LPC (only used for no PL) }
       inlag:integer;          { (i) pitch lag }
       iLBCdec_inst: piLBC_Dec_Inst_t
                           { (i/o) decoder instance }
   );
implementation
   procedure compCorr(
       cc:pareal;      { (o) cross correlation coefficient }
       gc:pareal;      { (o) gain }
       pm:pareal;
       buffer:PAreal;  { (i) signal buffer }
       lag:integer;    { (i) pitch lag }
       bLen:integer;       { (i) length of buffer }
       sRange:integer      { (i) correlation search length }
   );
   var
       i:integer;
       ftmp1, ftmp2, ftmp3:real;
   begin

       { Guard against getting outside buffer }
       if ((bLen-sRange-lag)<0) then
       begin
           sRange:=bLen-lag;
       end;

       ftmp1 := 0.0;
       ftmp2 := 0.0;
       ftmp3 := 0.0;
       for i:=0 to sRange-1 do
       begin
           ftmp1 :=ftmp1+ buffer[bLen-sRange+i]     * buffer[bLen-sRange+i-lag];
           ftmp2 :=ftmp2+ buffer[bLen-sRange+i-lag] * buffer[bLen-sRange+i-lag];
           ftmp3 :=ftmp3+ buffer[bLen-sRange+i]     * buffer[bLen-sRange+i];
       end;

       if (ftmp2 > 0.0) then
       begin
           cc[0] := ftmp1*ftmp1/ftmp2;
           gc[0] := abs(ftmp1/ftmp2);
           pm[0] :=abs(ftmp1)/(sqrt(ftmp2)*sqrt(ftmp3));
       end
       else
       begin
           cc[0] := 0.0;
           gc[0] := 0.0;
           pm[0] := 0.0;
       end;
   end;

   {----------------------------------------------------------------*
    *  Packet loss concealment routine. Conceals a residual signal
    *  and LP parameters. If no packet loss, update state.
    *---------------------------------------------------------------}

   procedure doThePLC(
       PLCresidual:pareal; { (o) concealed residual }
       PLClpc:pareal;      { (o) concealed LP parameters }
       PLI:integer;        { (i) packet loss indicator
                                  0 - no PL, 1 = PL }
       decresidual:pareal; { (i) decoded residual }
       lpc:pareal;         { (i) decoded LPC (only used for no PL) }
       inlag:integer;          { (i) pitch lag }
       iLBCdec_inst: piLBC_Dec_Inst_t
                           { (i/o) decoder instance }
   );
   var
     lag, randlag:integer;
     gain, maxcc:real;
     use_gain:real;
     gain_comp, maxcc_comp, per, max_per:real;
     i, pick, use_lag:integer;
     ftmp, pitchfact, energy:real;
     randvec:array [0..BLOCKL_MAX-1] of real;
   begin
			lag:=20;

       { Packet Loss }

       if (PLI = 1) then
       begin

           iLBCdec_inst^.consPLICount := iLBCdec_inst^.consPLICount+1;

           { if previous frame not lost,
              determine pitch pred. gain }

           if (iLBCdec_inst^.prevPLI <> 1) then
           begin

               { Search around the previous lag to find the
                  best pitch period }
               lag:=inlag-3;
               compCorr(@maxcc, @gain, @max_per,@iLBCdec_inst^.prevResidual,lag, iLBCdec_inst^.blockl, 60);
               for i:=inlag-2 to inlag+3 do
               begin
                   compCorr(@maxcc_comp, @gain_comp, @per,@iLBCdec_inst^.prevResidual,i, iLBCdec_inst^.blockl, 60);

                   if (maxcc_comp>maxcc) then
                   begin
                       maxcc:=maxcc_comp;
                       gain:=gain_comp;
                       lag:=i;
                       max_per:=per;
                   end;
               end;
          end

           { previous frame lost, use recorded lag and periodicity }

           else 
           begin
               lag:=iLBCdec_inst^.prevLag;
               max_per:=iLBCdec_inst^.per;
           end;

           { downscaling }

           use_gain:=1.0;
           if (iLBCdec_inst^.consPLICount*iLBCdec_inst^.blockl>320) then
               use_gain:=0.9
           else 
           if (iLBCdec_inst^.consPLICount*iLBCdec_inst^.blockl>2*320) then
               use_gain:=0.7
           else 
           if (iLBCdec_inst^.consPLICount*iLBCdec_inst^.blockl>3*320) then
               use_gain:=0.5
           else 
           if (iLBCdec_inst^.consPLICount*iLBCdec_inst^.blockl>4*320) then
               use_gain:=0.0;

           { mix noise and pitch repeatition }
           ftmp:=sqrt(max_per);
           if (ftmp> 0.7) then
               pitchfact:=1.0
           else 
           if (ftmp>0.4) then
               pitchfact:=(ftmp-0.4)/(0.7-0.4)
           else
               pitchfact:=0.0;


           { avoid repetition of same pitch cycle }
           use_lag:=lag;
           if (lag<80) then
           begin
               use_lag:=2*lag;
           end;

           { compute concealed residual }

           energy := 0.0;
           for i:=0 to iLBCdec_inst^.blockl-1 do
           begin
               { noise component }
               iLBCdec_inst^.seed:=(iLBCdec_inst^.seed*69069+1) and ($80000000-1);
               randlag := 50 + (iLBCdec_inst^.seed) mod 70;
               pick := i - randlag;

               if (pick < 0) then
               begin
                   randvec[i] := iLBCdec_inst^.prevResidual[iLBCdec_inst^.blockl+pick];
               end
               else
               begin
                   randvec[i] :=  randvec[pick];
               end;

               { pitch repeatition component }
               pick := i - use_lag;

               if (pick < 0) then
               begin
                   PLCresidual[i] :=iLBCdec_inst^.prevResidual[iLBCdec_inst^.blockl+pick];
               end
               else
               begin
                   PLCresidual[i] := PLCresidual[pick];
               end;

               { mix random and periodicity component }

               if (i<80) then
                   PLCresidual[i] := use_gain*(pitchfact * PLCresidual[i] + (1.0 - pitchfact) * randvec[i])
               else
               if (i<160) then
                   PLCresidual[i] := 0.95*use_gain*(pitchfact *PLCresidual[i] +(1.0 - pitchfact) * randvec[i])
               else
                   PLCresidual[i] := 0.9*use_gain*(pitchfact *PLCresidual[i] +(1.0 - pitchfact) * randvec[i]);

               energy :=energy + PLCresidual[i] * PLCresidual[i];
           end;

           { less than 30 dB, use only noise }






           if (sqrt(energy/iLBCdec_inst^.blockl) < 30.0) then
           begin
               gain:=0.0;
               for i:=0 to iLBCdec_inst^.blockl-1 do
               begin
                   PLCresidual[i] := randvec[i];
               end;
           end;

           { use old LPC }

           move(iLBCdec_inst^.prevLpc[0],PLClpc[0],(LPC_FILTERORDER+1)*sizeof(real));
       end

       { no packet loss, copy input }

       else
       begin
           move( decresidual[0],PLCresidual[0],iLBCdec_inst^.blockl*sizeof(real));
           move( lpc[0],PLClpc[0], (LPC_FILTERORDER+1)*sizeof(real));
           iLBCdec_inst^.consPLICount := 0;
       end;

       { update state }

       if (PLI=1) then
       begin
           iLBCdec_inst^.prevLag := lag;
           iLBCdec_inst^.per     := max_per;
       end;

       iLBCdec_inst^.prevPLI := PLI;
       move ( PLClpc[0],iLBCdec_inst^.prevLpc[0], (LPC_FILTERORDER+1)*sizeof(real));
       move ( PLCresidual[0],iLBCdec_inst^.prevResidual[0],iLBCdec_inst^.blockl*sizeof(real));
   end;
end.
