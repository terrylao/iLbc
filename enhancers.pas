{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit enhancers;

{$MODE Delphi}

interface
uses iLBC_define,constants,filter,C2Delphi_header;

   {----------------------------------------------------------------*
    * Find index in array such that the array element with said
    * index is the element of said array closest to "value"
    * according to the squared-error criterion
    *---------------------------------------------------------------}
   procedure NearestNeighbor(
       index:pInteger;   { (o) index of array element closest
                              to value }
       narray:PAreal;   { (i) data array }
       value:real;{ (i) value }
       arlength:integer{ (i) dimension of data array }
   );
   procedure mycorr1(
       corr:pareal;    { (o) correlation of seq1 and seq2 }
       seq1:pareal;    { (i) first sequence }
       dim1:integer;           { (i) dimension first seq1 }
       seq2:pareal;  { (i) second sequence }
       dim2:integer        { (i) dimension seq2 }
   );
   procedure enh_upsample(
       useq1:pareal;   { (o) upsampled output sequence }
       seq1:pareal;{ (i) unupsampled sequence }
       dim1:integer;       { (i) dimension seq1 }
       hfl:integer         { (i) polyphase filter length=2*hfl+1 }
   );
   procedure refiner(
       seg:pareal;         { (o) segment array }
       updStartPos:pareal; { (o) updated start point }
       idata:pareal;       { (i) original data buffer }
       idatal:integer;         { (i) dimension of idata }
       centerStartPos:integer; { (i) beginning center segment }
       estSegPos:real;{ (i) estimated beginning other segment }
       period:real    { (i) estimated pitch period }
   );
   procedure smath(
       odata:PAreal;   { (o) smoothed output }
       sseq:PAreal;{ (i) said second sequence of waveforms }
       hl:integer;         { (i) 2*hl+1 is sseq dimension }
       alpha0:real{ (i) max smoothing energy fraction }
   );
   procedure getsseq(
       sseq:PAreal;   { (o) the pitch-synchronous sequence }
       idata:PAreal;       { (i) original data }
       idatal:integer;         { (i) dimension of data }
       centerStartPos:integer; { (i) where current block starts }
       period:pareal;      { (i) rough-pitch-period array }
       plocs:pareal;       { (i) where periods of period array
                                  are taken }
       periodl:integer;    { (i) dimension period array }
       hl:integer              { (i) 2*hl+1 is the number of sequences }
   );
   procedure enhancer(
       odata:pareal;       { (o) smoothed block, dimension blockl }
       idata:pareal;       { (i) data buffer used for enhancing }
       idatal:integer;         { (i) dimension idata }
       centerStartPos:integer; { (i) first sample current block
                                  within idata }
       alpha0:real;       { (i) max correction-energy-fraction
                                 (in [0,1]) }
       period:pareal;      { (i) pitch period array }
       plocs:pareal;       { (i) locations where period array
                                  values valid }
       periodl:integer         { (i) dimension of period and plocs }
   );
   function xCorrCoef(
       target:PAreal;      { (i) first array }
       regressor:PAreal;   { (i) second array }
       subl:integer        { (i) dimension arrays }
   ):real;
   function enhancerInterface(
       nout:PAreal;                     { (o) enhanced signal }
       nin:PAreal;                      { (i) unenhanced signal }
       iLBCdec_inst:piLBC_Dec_Inst_t   { (i) buffers etc }
   ):integer;
implementation
   procedure NearestNeighbor(
       index:pInteger;   { (o) index of array element closest
                              to value }
       narray:PAreal;   { (i) data array }
       value:real;{ (i) value }
       arlength:integer{ (i) dimension of data array }
   );
   var
	   i:integer;
	   bestcrit,crit:real;
   begin
       crit:=narray[0]-value;
       bestcrit:=crit*crit;
       index^:=0;
       for i:=1 to arlength-1 do
       begin
           crit:=narray[i]-value;
           crit:=crit*crit;

           if (crit<bestcrit) then
           begin
               bestcrit:=crit;
               index^:=i;
           end;
       end;
   end;

   {----------------------------------------------------------------*
    * compute cross correlation between sequences
    *---------------------------------------------------------------}

   procedure mycorr1(
       corr:pareal;    { (o) correlation of seq1 and seq2 }
       seq1:pareal;    { (i) first sequence }
       dim1:integer;           { (i) dimension first seq1 }
       seq2:pareal;  { (i) second sequence }
       dim2:integer        { (i) dimension seq2 }
   );
   var
     i,j:integer;
   begin
       for i:=0 to dim1-dim2 do
       begin
           corr[i]:=0.0;
           for j:=0 to dim2-1 do
           begin
               corr[i] := corr[i] + seq1[i+j] * seq2[j];
           end;
       end;
   end;

   {----------------------------------------------------------------*
    * upsample finite array assuming zeros outside bounds
    *---------------------------------------------------------------}

   procedure enh_upsample(
       useq1:pareal;   { (o) upsampled output sequence }
       seq1:pareal;{ (i) unupsampled sequence }
       dim1:integer;       { (i) dimension seq1 }
       hfl:integer         { (i) polyphase filter length=2*hfl+1 }
   );
   var
     pu,ps:preal;
     i,j,k,q,filterlength,hfl2:integer;
     polyp:array [0..ENH_UPS0-1] of preal; { pointers to
                                      polyphase columns }
     pp:preal;
   begin
       { define pointers for filter }
       filterlength:=2*hfl+1;
       if ( filterlength > dim1 ) then
       begin
           hfl2:=(dim1 div 2);
           for j:=0 to ENH_UPS0-1 do
           begin
               polyp[j]:=@polyphaserTbl[j*filterlength+hfl-hfl2];
           end;
           hfl:=hfl2;
           filterlength:=2*hfl+1;
       end
       else
       begin
           for j:=0 to ENH_UPS0-1 do
           begin
               polyp[j]:=@polyphaserTbl[j*filterlength];
           end;
       end;

       { filtering: filter overhangs left side of sequence }

       pu:=@useq1[0];
       for i:=hfl to filterlength-1 do
       begin
           for j:=0 to ENH_UPS0-1 do
           begin
               pu^:=0.0;
               pp := polyp[j];
               ps := @seq1[i];
               for k:=0 to i do
               begin
                   pu^ := pu^ + ps^ * pp^;
                   dec(ps);
                   inc(pp);
               end;
               inc(pu);
           end;
       end;

       { filtering: simple convolution=inner products }

       for i:=filterlength to dim1-1 do
       begin
           for j:=0 to ENH_UPS0-1 do
           begin
               pu^:=0.0;
               pp := polyp[j];
               ps := @seq1[i];
               for k:=0 to filterlength-1 do
               begin
                   pu^ :=pu^ + ps^ * pp^;
                   dec(ps);
                   inc(pp);
               end;
               inc(pu);
           end;
       end;

       { filtering: filter overhangs right side of sequence }

       for q:=1 to hfl do
       begin
           for j:=0 to ENH_UPS0-1 do
           begin
               pu^:=0.0;
               inc(polyp[j],q);
               pp := polyp[j];
               dec(polyp[j],q);
               ps := @seq1[dim1-1];
               for k:=0 to filterlength-q-1 do
               begin
                   pu^ := pu^ + ps^ * pp^;
                   dec(ps);
                   inc(pp);
               end;
               inc(pu);
           end;
       end;
   end;


   {----------------------------------------------------------------*
    * find segment starting near idata+estSegPos that has highest
    * correlation with idata+centerStartPos through
    * idata+centerStartPos+ENH_BLOCKL-1 segment is found at a
    * resolution of ENH_UPSO times the original of the original
    * sampling rate
    *---------------------------------------------------------------}

   procedure refiner(
       seg:pareal;         { (o) segment array }
       updStartPos:pareal; { (o) updated start point }
       idata:pareal;       { (i) original data buffer }
       idatal:integer;         { (i) dimension of idata }
       centerStartPos:integer; { (i) beginning center segment }
       estSegPos:real;{ (i) estimated beginning other segment }
       period:real    { (i) estimated pitch period }
   );
   var
     estSegPosRounded,searchSegStartPos,searchSegEndPos,corrdim:integer;
     tloc,tloc2,i,st,en,fraction:integer;
     vect:array [0..ENH_VECTL-1] of real;
     corrVec:array [0..ENH_CORRDIM-1] of real;
     maxv:real;
     corrVecUps:array [0..ENH_CORRDIM*ENH_UPS0-1] of real;
   begin
       { defining array bounds }

       estSegPosRounded:=Trunc((estSegPos - 0.5));

       searchSegStartPos:=estSegPosRounded-ENH_SLOP;

       if (searchSegStartPos<0) then
       begin
           searchSegStartPos:=0;
       end;
       searchSegEndPos:=estSegPosRounded+ENH_SLOP;

       if (searchSegEndPos+ENH_BLOCKL >= idatal) then
       begin
           searchSegEndPos:=idatal-ENH_BLOCKL-1;
       end;
       corrdim:=searchSegEndPos-searchSegStartPos+1;

       { compute upsampled correlation (corr33) and find
          location of max }

       mycorr1(@corrVec,@idata[searchSegStartPos],
           corrdim+ENH_BLOCKL-1,@idata[centerStartPos],ENH_BLOCKL);
       enh_upsample(@corrVecUps[0],@corrVec[0],corrdim,ENH_FL0);
       tloc:=0;
       maxv:=corrVecUps[0];
       for i:=1 to ENH_UPS0*corrdim-1 do
       begin
           if (corrVecUps[i]>maxv) then
           begin
               tloc:=i;
               maxv:=corrVecUps[i];
           end;
       end;

       { make vector can be upsampled without ever running outside
          bounds }

       updStartPos[0]:= searchSegStartPos + tloc/ENH_UPS0+1.0;
       tloc2:=trunc(tloc/ENH_UPS0);

       if (tloc>tloc2*ENH_UPS0) then
       begin
           inc(tloc2);
       end;
       st:=searchSegStartPos+tloc2-ENH_FL0;

       if (st<0) then
       begin
           fillchar(vect,-st*sizeof(real),0);
           move(idata[0],vect[-st], (ENH_VECTL+st)*sizeof(real));
       end
       else
       begin
           en:=st+ENH_VECTL;

           if (en>idatal) then
           begin
               move( idata[st],vect[0],(ENH_VECTL-(en-idatal))*sizeof(real));
               fillchar(vect[ENH_VECTL-(en-idatal)],(en-idatal)*sizeof(real), 0);
           end
           else
           begin
               move( idata[st],vect[0], ENH_VECTL*sizeof(real));
           end;
       end;
       fraction:=tloc2*ENH_UPS0-tloc;

       { compute the segment (this is actually a convolution) }

       mycorr1(seg,@vect[0],ENH_VECTL,@polyphaserTbl[(2*ENH_FL0+1)*fraction],
           2*ENH_FL0+1);
   end;

   {----------------------------------------------------------------*
    * find the smoothed output data
    *---------------------------------------------------------------}

   procedure smath(
       odata:PAreal;   { (o) smoothed output }
       sseq:PAreal;{ (i) said second sequence of waveforms }
       hl:integer;         { (i) 2*hl+1 is sseq dimension }
       alpha0:real{ (i) max smoothing energy fraction }
   );
   var
     i,k:integer;
     psseq:PAreal;
     w00,w10,w11,A,B,C,err,errs:real;
     surround:array [0..BLOCKL_MAX-1] of real; { shape contributed by other than
                                    current }
     wt:array [0..2*ENH_HL+1 -1] of real;       { waveform weighting to get
                                    surround shape }
     denom:real;
   begin
       { create shape of contribution from all waveforms except the
          current one }

       for i:=1 to 2*hl+1 do
       begin
           wt[i-1] := 0.5*(1 - cos(2*PI*i/(2*hl+2)));
       end;
       wt[hl]:=0.0; { for clarity, not used }
       for i:=0 to ENH_BLOCKL-1 do
       begin
           surround[i]:=sseq[i]*wt[0];
       end;

       for k:=1 to hl-1 do
       begin
           psseq:=@sseq[k*ENH_BLOCKL];
           for i:=0 to ENH_BLOCKL-1 do
           begin
               surround[i]:=surround[i]+psseq[i]*wt[k];
           end;
       end;
       for k:=hl+1 to 2*hl do
       begin
           psseq:=@sseq[k*ENH_BLOCKL];
           for i:=0 to ENH_BLOCKL-1 do
           begin
               surround[i]:=surround[i]+psseq[i]*wt[k];
           end;
       end;

       { compute some inner products }

       w00 :=0.0;
       w10 :=0.0;
       w11 :=0.0;
       psseq:=@sseq[hl*ENH_BLOCKL]; { current block  }
       for i:=0 to ENH_BLOCKL-1 do
       begin
           w00:=w00+psseq[i]*psseq[i];
           w11:=w11+surround[i]*surround[i];
           w10:=w10+surround[i]*psseq[i];
       end;

       if (abs(w11) < 1.0) then
       begin
           w11:=1.0;
       end;
       C := sqrt( w00/w11);

       { first try enhancement without power-constraint }

       errs:=0.0;
       psseq:=@sseq[hl*ENH_BLOCKL];
       for i:=0 to ENH_BLOCKL-1 do
       begin
           odata[i]:=C*surround[i];
           err:=psseq[i]-odata[i];
           errs:=errs+err*err;
       end;

       { if constraint violated by first try, add constraint }

       if (errs > alpha0 * w00) then
       begin
           if ( w00 < 1) then
           begin
               w00:=1;
           end;
           denom := (w11*w00-w10*w10)/(w00*w00);

           if (denom > 0.0001) then { eliminates numerical problems for if smooth }
           begin
               A := sqrt( (alpha0- alpha0*alpha0/4)/denom);
               B := -alpha0/2 - A * w10/w00;
               B := B+1;
           end
           else
           begin { essentially no difference between cycles;
                     smoothing not needed }
               A:= 0.0;
               B:= 1.0;
           end;

           { create smoothed sequence }

           psseq:=@sseq[hl*ENH_BLOCKL];
           for i:=0 to ENH_BLOCKL-1 do
           begin
               odata[i]:=A*surround[i]+B*psseq[i];
           end;
       end;
   end;

   {----------------------------------------------------------------*
    * get the pitch-synchronous sample sequence
    *---------------------------------------------------------------}

   procedure getsseq(
       sseq:PAreal;   { (o) the pitch-synchronous sequence }
       idata:PAreal;       { (i) original data }
       idatal:integer;         { (i) dimension of data }
       centerStartPos:integer; { (i) where current block starts }
       period:pareal;      { (i) rough-pitch-period array }
       plocs:pareal;       { (i) where periods of period array
                                  are taken }
       periodl:integer;    { (i) dimension period array }
       hl:integer              { (i) 2*hl+1 is the number of sequences }
   );
   var
       i,centerEndPos,q:integer;
       blockStartPos:array [0..2*ENH_HL+1-1] of real;
       lagBlock:array [0..2*ENH_HL+1-1] of integer;
       plocs2:array [0..ENH_PLOCSL-1] of real;
       psseq:^real;
   begin
       centerEndPos:=centerStartPos+ENH_BLOCKL-1;

       { present }

       NearestNeighbor(@lagBlock[hl],plocs,0.5*(centerStartPos+centerEndPos),periodl);

       blockStartPos[hl]:=centerStartPos;

       psseq:=@sseq[ENH_BLOCKL*hl];
       move(idata[centerStartPos],psseq^,  ENH_BLOCKL*sizeof(real));

       { past }

       for q:=hl-1 downto 0 do
       begin
           blockStartPos[q]:=blockStartPos[q+1]-period[lagBlock[q+1]];
           NearestNeighbor(@lagBlock[q],plocs,
               blockStartPos[q]+ENH_BLOCKL_HALF-period[lagBlock[q+1]], periodl);


           if (blockStartPos[q]-ENH_OVERHANG>=0) then
           begin
               refiner(@sseq[q*ENH_BLOCKL], @blockStartPos[q], idata,
                   idatal, centerStartPos, blockStartPos[q],
                   period[lagBlock[q+1]]);
           end
           else
           begin
               psseq:=@sseq[q*ENH_BLOCKL];
               fillchar(psseq^, ENH_BLOCKL*sizeof(real), 0);
           end;
       end;

       { future }

       for i:=0 to periodl-1 do
       begin
           plocs2[i]:=plocs[i]-period[i];
       end;
       for q:=hl+1 to 2*hl do
       begin
           NearestNeighbor(@lagBlock[q],@plocs2,
               blockStartPos[q-1]+ENH_BLOCKL_HALF,periodl);

           blockStartPos[q]:=blockStartPos[q-1]+period[lagBlock[q]];
           if (blockStartPos[q]+ENH_BLOCKL+ENH_OVERHANG<idatal) then
           begin
               refiner(@sseq[ENH_BLOCKL*q], @blockStartPos[q], idata,
                   idatal, centerStartPos, blockStartPos[q],
                   period[lagBlock[q]]);
           end
           else
           begin
               psseq:=@sseq[q*ENH_BLOCKL];
               fillchar(psseq^, ENH_BLOCKL*sizeof(real), 0);
           end;
       end;
       
   end;

   {----------------------------------------------------------------*
    * perform enhancement on idata+centerStartPos through
    * idata+centerStartPos+ENH_BLOCKL-1
    *---------------------------------------------------------------}

   procedure enhancer(
       odata:pareal;       { (o) smoothed block, dimension blockl }
       idata:pareal;       { (i) data buffer used for enhancing }
       idatal:integer;         { (i) dimension idata }
       centerStartPos:integer; { (i) first sample current block
                                  within idata }
       alpha0:real;       { (i) max correction-energy-fraction
                                 (in [0,1]) }
       period:pareal;      { (i) pitch period array }
       plocs:pareal;       { (i) locations where period array
                                  values valid }
       periodl:integer         { (i) dimension of period and plocs }
   );
   var
   	sseq:array [0..(2*ENH_HL+1)*ENH_BLOCKL-1] of real;
   begin
       { get said second sequence of segments }

       getsseq(@sseq,idata,idatal,centerStartPos,period,
           plocs,periodl,ENH_HL);

       { compute the smoothed output from said second sequence }

       smath(odata,@sseq,ENH_HL,alpha0);
  end;

   {----------------------------------------------------------------*
    * cross correlation
    *---------------------------------------------------------------}

   function xCorrCoef(
       target:PAreal;      { (i) first array }
       regressor:PAreal;   { (i) second array }
       subl:integer        { (i) dimension arrays }
   ):real;
   var
       i:integer;
       ftmp1, ftmp2:real;
   begin
       ftmp1 := 0.0;
       ftmp2 := 0.0;
       for i:=0 to subl-1 do
       begin
           ftmp1 :=ftmp1 + target[i]*regressor[i];
           ftmp2 :=ftmp2 + regressor[i]*regressor[i];
       end;

       if (ftmp1 > 0.0) then
       begin
           result:=(ftmp1*ftmp1/ftmp2);
       end
       else
       begin
           result:=0.0;
       end;
   end;

   {----------------------------------------------------------------*
    * interface for enhancer
    *---------------------------------------------------------------}
   function enhancerInterface(
       nout:PAreal;                     { (o) enhanced signal }
       nin:PAreal;                      { (i) unenhanced signal }
       iLBCdec_inst:piLBC_Dec_Inst_t   { (i) buffers etc }
   ):integer;
   var
       enh_buf,enh_period:PAreal;
       iblock, isample:integer;
       lag, ilag, i, ioffset:integer;
       cc, maxcc:real;
       ftmp1, ftmp2:real;
       inPtr, enh_bufPtr1, enh_bufPtr2:^real;
       plc_pred:array [0..ENH_BLOCKL-1] of real;

       lpState:array[0..5] of real;
       downsampled:array [0..((ENH_NBLOCKS*ENH_BLOCKL+120) div 2) - 1] of real;
       inLen:integer;
       start, plc_blockl, inlag:integer;
   begin
   		lag:=0;
   		inLen:=ENH_NBLOCKS*ENH_BLOCKL+120;
       enh_buf:=@iLBCdec_inst^.enh_buf;
       enh_period:=@iLBCdec_inst^.enh_period;
       move(enh_buf[iLBCdec_inst^.blockl],enh_buf[0], (ENH_BUFL-iLBCdec_inst^.blockl)*sizeof(real));

       move(nin[0],enh_buf[ENH_BUFL-iLBCdec_inst^.blockl], iLBCdec_inst^.blockl*sizeof(real));

       if (iLBCdec_inst^.mode=30) then
           plc_blockl:=ENH_BLOCKL
       else
           plc_blockl:=40;

       { when 20 ms frame, move processing one block }
       ioffset:=0;
       if (iLBCdec_inst^.mode=20) then
       	ioffset:=1;

       i:=3-ioffset;
       move(enh_period[i], enh_period[0], (ENH_NBLOCKS_TOT-i)*sizeof(real));

       { Set state information to the 6 samples right before
          the samples to be downsampled. }

       move(enh_buf[(ENH_NBLOCKS_EXTRA+ioffset)*ENH_BLOCKL-126],lpState,6*sizeof(real));

       { Down sample a factor 2 to save computations }

       DownSample(@enh_buf[(ENH_NBLOCKS_EXTRA+ioffset)*ENH_BLOCKL-120],
                   @lpFilt_coefsTbl, inLen-ioffset*ENH_BLOCKL,
                   @lpState, @downsampled);

       { Estimate the pitch in the down sampled domain. }
       for iblock := 0 to ENH_NBLOCKS-ioffset-1 do
			 begin
           lag := 10;
           maxcc := xCorrCoef(@downsampled[60+iblock*ENH_BLOCKL_HALF], @downsampled[60+iblock*ENH_BLOCKL_HALF-lag], ENH_BLOCKL_HALF);
           for ilag:=11 to 59 do
           begin
               cc := xCorrCoef(@downsampled[60+iblock*
                   ENH_BLOCKL_HALF], @downsampled[60+iblock*
                   ENH_BLOCKL_HALF-ilag], ENH_BLOCKL_HALF);

               if (cc > maxcc) then
               begin
                   maxcc := cc;
                   lag := ilag;
               end;
           end;

           { Store the estimated lag in the non-downsampled domain }
           enh_period[iblock+ENH_NBLOCKS_EXTRA+ioffset] := lag*2;
       end;


       { PLC was performed on the previous packet }
       if (iLBCdec_inst^.prev_enh_pl=1) then
       begin
           inlag:=trunc(enh_period[ENH_NBLOCKS_EXTRA+ioffset]);

           lag := inlag-1;
           maxcc := xCorrCoef(nin, @nin[lag], plc_blockl);
           for ilag:=inlag to inlag+1 do
           begin
               cc := xCorrCoef(nin, @nin[ilag], plc_blockl);
               if (cc > maxcc) then
               begin
                   maxcc := cc;
                   lag := ilag;
               end;
           end;

           enh_period[ENH_NBLOCKS_EXTRA+ioffset-1]:=lag;

           { compute new concealed residual for the old lookahead,
              mix the forward PLC with a backward PLC from
              the new frame }

           inPtr:=@nin[lag-1];

           enh_bufPtr1:=@plc_pred[plc_blockl-1];

           if (lag>plc_blockl) then
           begin
               start:=plc_blockl;
           end
           else
           begin
               start:=lag;
           end;

           for isample := start downto 1 do
           begin
               enh_bufPtr1^ := inPtr^;
               dec(enh_bufPtr1);
               dec(inPtr);
           end;

           enh_bufPtr2:=@enh_buf[ENH_BUFL-1-iLBCdec_inst^.blockl];
           for isample := (plc_blockl-1-lag) downto 0 do
           begin
               enh_bufPtr1^ := enh_bufPtr2^;
               dec(enh_bufPtr1);
               dec(enh_bufPtr2);
           end;

           { limit energy change }
           ftmp2:=0.0;
           ftmp1:=0.0;
           for i:=0 to plc_blockl-1 do
           begin
               ftmp2:=ftmp2+enh_buf[ENH_BUFL-1-iLBCdec_inst^.blockl-i]*enh_buf[ENH_BUFL-1-iLBCdec_inst^.blockl-i];
               ftmp1:=ftmp1+plc_pred[i]*plc_pred[i];
           end;
           ftmp1:=sqrt(ftmp1/plc_blockl);
           ftmp2:=sqrt(ftmp2/plc_blockl);
           if (ftmp1>2.0*ftmp2) and (ftmp1>0.0) then
           begin
               for i:=0 to plc_blockl-9 do
               begin
                   plc_pred[i]:=plc_pred[i]*2.0*ftmp2/ftmp1;
               end;
               for i:=plc_blockl-10 to plc_blockl-1 do
               begin
                   plc_pred[i]:=plc_pred[i]*(i-plc_blockl+10)*(1.0-2.0*ftmp2/ftmp1)/(10)+2.0*ftmp2/ftmp1;
               end;
           end;

           enh_bufPtr1:=@enh_buf[ENH_BUFL-1-iLBCdec_inst^.blockl];
           for i:=0 to plc_blockl-1 do
           begin
               ftmp1 := (i+1) / (plc_blockl+1);
               enh_bufPtr1^ :=enh_bufPtr1^ * ftmp1;
               enh_bufPtr1^ :=enh_bufPtr1^ + (1.0-ftmp1)*plc_pred[plc_blockl-1-i];
               dec(enh_bufPtr1);
           end;
       end;

       if (iLBCdec_inst^.mode=20) then
       begin
           { Enhancer with 40 samples delay }
           for iblock := 0 to 1 do
           begin
               enhancer(@nout[iblock*ENH_BLOCKL], enh_buf,
                   ENH_BUFL, (5+iblock)*ENH_BLOCKL+40,
                   ENH_ALPHA0, enh_period, @enh_plocsTbl,
                       ENH_NBLOCKS_TOT);
           end;
       end
       else 
       if (iLBCdec_inst^.mode=30) then
       begin
           { Enhancer with 80 samples delay }
           for iblock := 0 to 2 do
           begin
               enhancer(@nout[iblock*ENH_BLOCKL], enh_buf,
                   ENH_BUFL, (4+iblock)*ENH_BLOCKL,
                   ENH_ALPHA0, enh_period, @enh_plocsTbl,
                       ENH_NBLOCKS_TOT);
           end;
       end;

       result:=lag*2;
   end;
end.
