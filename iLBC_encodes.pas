{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit iLBC_encodes;

{$MODE Delphi}

interface
uses
iLBC_define,
LPCencodes,
FrameClassifys,
StateSearchWs,
StateConstructWs,
helpfun,
constants,
packing,
iCBSearchs,
iCBConstructs,
hpInputs,
anaFilters,
syntFilters,C2Delphi_header;
   function initEncode(                   { (o) Number of bytes
                                              encoded }
       iLBCenc_inst:piLBC_Enc_Inst_t;  { (i/o) Encoder instance }
       mode:integer                    { (i) frame size mode }
   ):Smallint;
   procedure iLBC_encode(
       bytes:pchar;           { (o) encoded data bits iLBC }
       block:pareal;                   { (o) speech vector to
                                              encode }
       iLBCenc_inst:piLBC_Enc_Inst_t   { (i/o) the general encoder
                                              state }
   );
   {----------------------------------------------------------------*
    *  Initiation of encoder instance.
    *---------------------------------------------------------------}
implementation
   function initEncode(                   { (o) Number of bytes
                                              encoded }
       iLBCenc_inst:piLBC_Enc_Inst_t;  { (i/o) Encoder instance }
       mode:integer                    { (i) frame size mode }
   ):Smallint;
   begin
       iLBCenc_inst^.mode := mode;
       if (mode=30) then 
	     begin
           iLBCenc_inst^.blockl := BLOCKL_30MS;
           iLBCenc_inst^.nsub := NSUB_30MS;
           iLBCenc_inst^.nasub := NASUB_30MS;
           iLBCenc_inst^.lpc_n := LPC_N_30MS;
           iLBCenc_inst^.no_of_bytes := NO_OF_BYTES_30MS;
           iLBCenc_inst^.no_of_words := NO_OF_WORDS_30MS;
           iLBCenc_inst^.state_short_len:=STATE_SHORT_LEN_30MS;
           { ULP init }
           iLBCenc_inst^.ULP_inst:=@ULP_30msTbl;
       end
       else 
       if (mode=20) then 
	     begin
           iLBCenc_inst^.blockl := BLOCKL_20MS;
           iLBCenc_inst^.nsub := NSUB_20MS;
           iLBCenc_inst^.nasub := NASUB_20MS;
           iLBCenc_inst^.lpc_n := LPC_N_20MS;
           iLBCenc_inst^.no_of_bytes := NO_OF_BYTES_20MS;
           iLBCenc_inst^.no_of_words := NO_OF_WORDS_20MS;
           iLBCenc_inst^.state_short_len:=STATE_SHORT_LEN_20MS;
           { ULP init }
           iLBCenc_inst^.ULP_inst:=@ULP_20msTbl;
       end
       else 
       begin
           result:=2;
           exit;
       end;

       fillchar(iLBCenc_inst^.anaMem,LPC_FILTERORDER*sizeof(real), 0);
       move( lsfmeanTbl[0],iLBCenc_inst^.lsfold[0],LPC_FILTERORDER*sizeof(real));
       move( lsfmeanTbl[0],iLBCenc_inst^.lsfdeqold[0],LPC_FILTERORDER*sizeof(real));
       fillchar(iLBCenc_inst^.lpc_buffer[0],(LPC_LOOKBACK+BLOCKL_MAX)*sizeof(real), 0);
       fillchar(iLBCenc_inst^.hpimem[0], 4*sizeof(real), 0);
       result:=Smallint(iLBCenc_inst^.no_of_bytes);
   end;

   {----------------------------------------------------------------*
    *  main encoder function
    *---------------------------------------------------------------}

   procedure iLBC_encode(
       bytes:pchar;           { (o) encoded data bits iLBC }
       block:pareal;                   { (o) speech vector to
                                              encode }
       iLBCenc_inst:piLBC_Enc_Inst_t   { (i/o) the general encoder
                                              state }
   );
   var
       data:array [0..BLOCKL_MAX-1] of real;
       residual:array [0..BLOCKL_MAX-1] of real;
       reverseResidual:array [0..BLOCKL_MAX-1] of real;
       start, idxForMax:integer; 
       idxVec:array [0..STATE_LEN-1] of integer;
       reverseDecresidual:array [0..BLOCKL_MAX-1] of real;
       mem:array [0..CB_MEML-1] of real;
       n, k, meml_gotten, Nfor, Nback, i, pos:integer;
       gain_index:array [0..CB_NSTAGES*NASUB_MAX-1] of integer;
       extra_gain_index:array [0..CB_NSTAGES-1] of integer;
       cb_index:array [0..CB_NSTAGES*NASUB_MAX-1] of integer;
       extra_cb_index:array [0..CB_NSTAGES-1] of integer;
       lsf_i:array [0..LSF_NSPLIT*LPC_N_MAX-1] of integer;
       pbytes:pchar;
       diff, start_pos, state_first:integer;
       en1, en2:real;
       index, ulp, firstpart:integer;
       subcount, subframe:integer;
       weightState:array [0..LPC_FILTERORDER-1] of real;
       syntdenum:array [0..NSUB_MAX*(LPC_FILTERORDER+1)-1] of real;
       weightdenum:array [0..NSUB_MAX*(LPC_FILTERORDER+1)-1] of real;
       decresidual:array [0..BLOCKL_MAX-1] of real;
   begin
       { high pass filtering of input signal if such is not done
              prior to calling this function }

       hpInput(block, iLBCenc_inst^.blockl, @data[0], @iLBCenc_inst^.hpimem);

       { otherwise simply copy }

       {memcpy(data,block,iLBCenc_inst^.blockl*sizeof(real));}

       { LPC of hp filtered input data }

       LPCencode(@syntdenum[0], @weightdenum[0], @lsf_i[0], @data[0], iLBCenc_inst);


       { inverse filter to get residual }

       for n:=0 to iLBCenc_inst^.nsub-1 do
	     begin
           anaFilter(@data[n*SUBL], @syntdenum[n*(LPC_FILTERORDER+1)], SUBL, @residual[n*SUBL], @iLBCenc_inst^.anaMem);
       end;

       { find state location }

       start := FrameClassify(iLBCenc_inst, @residual);

       { check if state should be in first or last part of the
       two subframes }

       diff := STATE_LEN - iLBCenc_inst^.state_short_len;
       en1 := 0;
       index := (start-1)*SUBL;

       for i := 0 to iLBCenc_inst^.state_short_len-1 do
	     begin
           en1 :=en1 + residual[index+i]*residual[index+i];
       end;
       en2 := 0;
       index := (start-1)*SUBL+diff;
       for i := 0 to iLBCenc_inst^.state_short_len-1 do
	     begin
           en2 :=en2 + residual[index+i]*residual[index+i];
       end;


       if (en1 > en2) then 
	     begin
           state_first := 1;
           start_pos := (start-1)*SUBL;
       end 
       else 
       begin
           state_first := 0;
           start_pos := (start-1)*SUBL + diff;
       end;

       { scalar quantization of state }

       StateSearchW(iLBCenc_inst, @residual[start_pos],
           @syntdenum[(start-1)*(LPC_FILTERORDER+1)],
           @weightdenum[(start-1)*(LPC_FILTERORDER+1)], @idxForMax,
           @idxVec, iLBCenc_inst^.state_short_len, state_first);

       StateConstructW(idxForMax, @idxVec,
           @syntdenum[(start-1)*(LPC_FILTERORDER+1)],
           @decresidual[start_pos], iLBCenc_inst^.state_short_len);

       { predictive quantization in state }

       if (state_first<>0) then 
	     begin { put adaptive part in the end }

           { setup memory }

           fillchar(mem,(CB_MEML-iLBCenc_inst^.state_short_len)*sizeof(real), 0);
           move(decresidual[start_pos],mem[CB_MEML-iLBCenc_inst^.state_short_len],iLBCenc_inst^.state_short_len*sizeof(real));
           fillchar(weightState, LPC_FILTERORDER*sizeof(real), 0);

           { encode sub-frames }

           iCBSearch(iLBCenc_inst, @extra_cb_index, @extra_gain_index,
               @residual[start_pos+iLBCenc_inst^.state_short_len],
               @mem[CB_MEML-stMemLTbl],
               stMemLTbl, diff, CB_NSTAGES,
               @weightdenum[start*(LPC_FILTERORDER+1)],
               @weightState, 0);

           { construct decoded vector }

           iCBConstruct(
               @decresidual[start_pos+iLBCenc_inst^.state_short_len],
               @extra_cb_index, @extra_gain_index,
               @mem[CB_MEML-stMemLTbl],
               stMemLTbl, diff, CB_NSTAGES);

       end
       else 
       begin { put adaptive part in the beginning }

           { create reversed vectors for prediction }

           for k:=0 to diff-1 do
	         begin
               reverseResidual[k] := residual[(start+1)*SUBL-1-(k+iLBCenc_inst^.state_short_len)];
           end;

           { setup memory }

           meml_gotten := iLBCenc_inst^.state_short_len;
           k:=0;
           while k < meml_gotten do
	         begin
               mem[CB_MEML-1-k] := decresidual[start_pos + k];
               inc(k);
           end;
           fillchar(mem[0], (CB_MEML-k)*sizeof(real), 0);
           fillchar(weightState[0], LPC_FILTERORDER*sizeof(real), 0);

           { encode sub-frames }

           iCBSearch(iLBCenc_inst, @extra_cb_index, @extra_gain_index,
               @reverseResidual, @mem[CB_MEML-stMemLTbl], stMemLTbl,
               diff, CB_NSTAGES,
               @weightdenum[(start-1)*(LPC_FILTERORDER+1)],
               @weightState, 0);

           { construct decoded vector }

           iCBConstruct(@reverseDecresidual, @extra_cb_index,
               @extra_gain_index, @mem[CB_MEML-stMemLTbl], stMemLTbl,
               diff, CB_NSTAGES);

           { get decoded residual from reversed vector }

           for k:=0 to diff-1 do
	         begin
               decresidual[start_pos-1-k] := reverseDecresidual[k];
           end;
       end;

       { counter for predicted sub-frames }

       subcount:=0;

       { forward prediction of sub-frames }

       Nfor := iLBCenc_inst^.nsub-start-1;


       if ( Nfor > 0 ) then 
	     begin

           { setup memory }

           fillchar(mem, (CB_MEML-STATE_LEN)*sizeof(real), 0);
           move( decresidual[(start-1)*SUBL],mem[CB_MEML-STATE_LEN],STATE_LEN*sizeof(real));
           fillchar(weightState, LPC_FILTERORDER*sizeof(real), 0);

           { loop over sub-frames to encode }

           for subframe:=0 to Nfor-1 do
	         begin
               { encode sub-frame }

               iCBSearch(iLBCenc_inst, @cb_index[subcount*CB_NSTAGES],
                   @gain_index[subcount*CB_NSTAGES],
                   @residual[(start+1+subframe)*SUBL],
                   @mem[CB_MEML-memLfTbl[subcount]],
                   memLfTbl[subcount], SUBL, CB_NSTAGES,
                   @weightdenum[(start+1+subframe)*
                               (LPC_FILTERORDER+1)],
                   @weightState, subcount+1);

               { construct decoded vector }

               iCBConstruct(@decresidual[(start+1+subframe)*SUBL],
                   @cb_index[subcount*CB_NSTAGES],
                   @gain_index[subcount*CB_NSTAGES],
                   @mem[CB_MEML-memLfTbl[subcount]],
                   memLfTbl[subcount], SUBL, CB_NSTAGES);

               { update memory }

               move( mem[SUBL], mem[0],(CB_MEML-SUBL)*sizeof(real));
               move(decresidual[(start+1+subframe)*SUBL],mem[CB_MEML-SUBL],SUBL*sizeof(real));
               fillchar(weightState, LPC_FILTERORDER*sizeof(real), 0);
               inc(subcount);
           end;
       end;


       { backward prediction of sub-frames }

       Nback := start-1;


       if ( Nback > 0 ) then 
	     begin

           { create reverse order vectors }

           for n:=0 to Nback-1 do
	         begin
               for k:=0 to SUBL-1 do
	             begin
                   reverseResidual[n*SUBL+k] :=residual[(start-1)*SUBL-1-n*SUBL-k];
                   reverseDecresidual[n*SUBL+k] :=decresidual[(start-1)*SUBL-1-n*SUBL-k];
               end;
           end;

           { setup memory }

           meml_gotten := SUBL*(iLBCenc_inst^.nsub+1-start);


           if ( meml_gotten > CB_MEML ) then 
	         begin
               meml_gotten:=CB_MEML;
           end;
           k:=0;
           while k < meml_gotten do
	         begin
               mem[CB_MEML-1-k] := decresidual[(start-1)*SUBL + k];
               inc(k);
           end;
           fillchar(mem, (CB_MEML-k)*sizeof(real), 0);
           fillchar(weightState, LPC_FILTERORDER*sizeof(real), 0);

           { loop over sub-frames to encode }

           for subframe:=0 to Nback-1 do
	         begin
               { encode sub-frame }

               iCBSearch(iLBCenc_inst, @cb_index[subcount*CB_NSTAGES],
                   @gain_index[subcount*CB_NSTAGES],
                   @reverseResidual[subframe*SUBL],
                   @mem[CB_MEML-memLfTbl[subcount]],
                   memLfTbl[subcount], SUBL, CB_NSTAGES,
                   @weightdenum[(start-2-subframe)*
                               (LPC_FILTERORDER+1)],
                   @weightState, subcount+1);

               { construct decoded vector }

               iCBConstruct(@reverseDecresidual[subframe*SUBL],
                   @cb_index[subcount*CB_NSTAGES],
                   @gain_index[subcount*CB_NSTAGES],
                   @mem[CB_MEML-memLfTbl[subcount]],
                   memLfTbl[subcount], SUBL, CB_NSTAGES);

               { update memory }

               move( mem[SUBL],mem[0], (CB_MEML-SUBL)*sizeof(real));
               move(reverseDecresidual[subframe*SUBL],mem[CB_MEML-SUBL],SUBL*sizeof(real));
               fillchar(weightState, LPC_FILTERORDER*sizeof(real), 0);
               inc(subcount);
           end;

           { get decoded residual from reversed vector }

           for i:=0 to SUBL*Nback-1 do
	         begin
               decresidual[SUBL*Nback - i - 1] :=reverseDecresidual[i];
           end;
       end;
       { end encoding part }

       { adjust index }
       index_conv_enc(@cb_index);

       { pack bytes }

       pbytes:=bytes;
       pos:=0;

       { loop over the 3 ULP classes }

       for ulp:=0 to 2 do
	     begin
           { LSF }
           for k:=0 to LSF_NSPLIT*iLBCenc_inst^.lpc_n-1 do
	         begin
               packsplit(@lsf_i[k], @firstpart, @lsf_i[k],
                   iLBCenc_inst^.ULP_inst^.lsf_bits[k][ulp],
                   iLBCenc_inst^.ULP_inst^.lsf_bits[k][ulp]+
                   iLBCenc_inst^.ULP_inst^.lsf_bits[k][ulp+1]+
                   iLBCenc_inst^.ULP_inst^.lsf_bits[k][ulp+2]);
               dopack( @pbytes, firstpart,
                   iLBCenc_inst^.ULP_inst^.lsf_bits[k][ulp], @pos);
           end;

           { Start block info }

           packsplit(@start, @firstpart, @start,
               iLBCenc_inst^.ULP_inst^.start_bits[ulp],
               iLBCenc_inst^.ULP_inst^.start_bits[ulp]+
               iLBCenc_inst^.ULP_inst^.start_bits[ulp+1]+
               iLBCenc_inst^.ULP_inst^.start_bits[ulp+2]);
           dopack( @pbytes, firstpart,
               iLBCenc_inst^.ULP_inst^.start_bits[ulp], @pos);

           packsplit(@state_first, @firstpart, @state_first,
               iLBCenc_inst^.ULP_inst^.startfirst_bits[ulp],
               iLBCenc_inst^.ULP_inst^.startfirst_bits[ulp]+
               iLBCenc_inst^.ULP_inst^.startfirst_bits[ulp+1]+
               iLBCenc_inst^.ULP_inst^.startfirst_bits[ulp+2]);
           dopack( @pbytes, firstpart,
               iLBCenc_inst^.ULP_inst^.startfirst_bits[ulp], @pos);

           packsplit(@idxForMax, @firstpart, @idxForMax,
               iLBCenc_inst^.ULP_inst^.scale_bits[ulp],
               iLBCenc_inst^.ULP_inst^.scale_bits[ulp]+
               iLBCenc_inst^.ULP_inst^.scale_bits[ulp+1]+
               iLBCenc_inst^.ULP_inst^.scale_bits[ulp+2]);
           dopack( @pbytes, firstpart,
               iLBCenc_inst^.ULP_inst^.scale_bits[ulp], @pos);

           for k:=0 to iLBCenc_inst^.state_short_len-1 do
	         begin
               packsplit(@idxVec[k], @firstpart, @idxVec[k],
                   iLBCenc_inst^.ULP_inst^.state_bits[ulp],
                   iLBCenc_inst^.ULP_inst^.state_bits[ulp]+
                   iLBCenc_inst^.ULP_inst^.state_bits[ulp+1]+
                   iLBCenc_inst^.ULP_inst^.state_bits[ulp+2]);
               dopack( @pbytes, firstpart,
                   iLBCenc_inst^.ULP_inst^.state_bits[ulp], @pos);
           end;

           { 23/22 (20ms/30ms) sample block }

           for k:=0 to CB_NSTAGES-1 do
	         begin
               packsplit(@extra_cb_index[k], @firstpart,
                   @extra_cb_index[k],
                   iLBCenc_inst^.ULP_inst^.extra_cb_index[k][ulp],
                   iLBCenc_inst^.ULP_inst^.extra_cb_index[k][ulp]+
                   iLBCenc_inst^.ULP_inst^.extra_cb_index[k][ulp+1]+
                   iLBCenc_inst^.ULP_inst^.extra_cb_index[k][ulp+2]);
               dopack( @pbytes, firstpart,
                   iLBCenc_inst^.ULP_inst^.extra_cb_index[k][ulp],
                   @pos);
           end;

           for k:=0 to CB_NSTAGES-1 do
	         begin
               packsplit(@extra_gain_index[k], @firstpart,
                   @extra_gain_index[k],
                   iLBCenc_inst^.ULP_inst^.extra_cb_gain[k][ulp],
                   iLBCenc_inst^.ULP_inst^.extra_cb_gain[k][ulp]+
                   iLBCenc_inst^.ULP_inst^.extra_cb_gain[k][ulp+1]+
                   iLBCenc_inst^.ULP_inst^.extra_cb_gain[k][ulp+2]);
               dopack( @pbytes, firstpart,
                   iLBCenc_inst^.ULP_inst^.extra_cb_gain[k][ulp],
                   @pos);
           end;

           { The two/four (20ms/30ms) 40 sample sub-blocks }

           for i:=0 to iLBCenc_inst^.nasub-1 do
	         begin
               for k:=0 to CB_NSTAGES-1 do
	             begin
                   packsplit(@cb_index[i*CB_NSTAGES+k], @firstpart,
                       @cb_index[i*CB_NSTAGES+k],
                       iLBCenc_inst^.ULP_inst^.cb_index[i][k][ulp],
                       iLBCenc_inst^.ULP_inst^.cb_index[i][k][ulp]+
                       iLBCenc_inst^.ULP_inst^.cb_index[i][k][ulp+1]+
                       iLBCenc_inst^.ULP_inst^.cb_index[i][k][ulp+2]);
                   dopack( @pbytes, firstpart,
                       iLBCenc_inst^.ULP_inst^.cb_index[i][k][ulp],
                       @pos);
               end;
           end;

           for i:=0 to iLBCenc_inst^.nasub-1 do
	         begin
               for k:=0 to CB_NSTAGES-1 do
	             begin
                   packsplit(@gain_index[i*CB_NSTAGES+k], @firstpart,
                       @gain_index[i*CB_NSTAGES+k],
                       iLBCenc_inst^.ULP_inst^.cb_gain[i][k][ulp],
                       iLBCenc_inst^.ULP_inst^.cb_gain[i][k][ulp]+
                       iLBCenc_inst^.ULP_inst^.cb_gain[i][k][ulp+1]+
                       iLBCenc_inst^.ULP_inst^.cb_gain[i][k][ulp+2]);
                   dopack( @pbytes, firstpart,
                       iLBCenc_inst^.ULP_inst^.cb_gain[i][k][ulp],
                       @pos);
               end;
           end;
       end;
       { set the last bit to zero (otherwise the decoder
          will treat it as a lost frame) }
       dopack( @pbytes, 0, 1, @pos);
   end;
end.
