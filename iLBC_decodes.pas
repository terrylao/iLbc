unit iLBC_decodes;

{$MODE Delphi}

interface
uses
iLBC_define,
StateConstructWs,
LPCdecode,
iCBConstructs,
doCPLC,
helpfun,
constants,
packing,
enhancers,
hpOutputs,
syntFilters,C2Delphi_header;

   {----------------------------------------------------------------*
    *  Initiation of decoder instance.
    *---------------------------------------------------------------}
   function initDecode(                   { (o) Number of decoded
                                              samples }
       iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) Decoder instance }
       mode:integer;                       { (i) frame size mode }
       use_enhancer:integer                { (i) 1 to use enhancer
                                              0 to run without
                                                enhancer }
   ):Smallint;
   procedure Decode(
       iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) the decoder state
                                                structure }
       decresidual:pareal;             { (o) decoded residual frame }
       start:integer;                      { (i) location of start
                                              state }
       idxForMax:integer;                  { (i) codebook index for the
                                              maximum value }
       idxVec:painteger;                { (i) codebook indexes for the
                                              samples  in the start
                                              state }
       syntdenum:pareal;               { (i) the decoded synthesis
                                              filter coefficients }
       cb_index:painteger;                  { (i) the indexes for the
                                              adaptive codebook }
       gain_index:painteger;            { (i) the indexes for the
                                              corresponding gains }
       extra_cb_index:painteger;        { (i) the indexes for the
                                              adaptive codebook part
                                              of start state }
       extra_gain_index:painteger;          { (i) the indexes for the
                                              corresponding gains }
       state_first:integer                 { (i) 1 if non adaptive part
                                              of start state comes
                                              first 0 if that part
                                              comes last }
   );
   procedure iLBC_decode(
       decblock:pareal;            { (o) decoded signal block }
       bytes:pchar;           { (i) encoded signal bits }
       iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) the decoder state
                                                structure }
       mode:integer                    { (i) 0: bad packet, PLC,
                                              1: normal }
   );
implementation
   function initDecode(                   { (o) Number of decoded
                                              samples }
       iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) Decoder instance }
       mode:integer;                       { (i) frame size mode }
       use_enhancer:integer                { (i) 1 to use enhancer
                                              0 to run without
                                                enhancer }
   ):Smallint;
   var
   	i:integer;
   begin

       iLBCdec_inst^.mode := mode;

       if (mode=30) then
       begin
           iLBCdec_inst^.blockl := BLOCKL_30MS;
           iLBCdec_inst^.nsub := NSUB_30MS;
           iLBCdec_inst^.nasub := NASUB_30MS;
           iLBCdec_inst^.lpc_n := LPC_N_30MS;
           iLBCdec_inst^.no_of_bytes := NO_OF_BYTES_30MS;
           iLBCdec_inst^.no_of_words := NO_OF_WORDS_30MS;
           iLBCdec_inst^.state_short_len:=STATE_SHORT_LEN_30MS;
           { ULP init }
           iLBCdec_inst^.ULP_inst:=@ULP_30msTbl;
       end
       else 
       if (mode=20) then
       begin
           iLBCdec_inst^.blockl := BLOCKL_20MS;
           iLBCdec_inst^.nsub := NSUB_20MS;
           iLBCdec_inst^.nasub := NASUB_20MS;
           iLBCdec_inst^.lpc_n := LPC_N_20MS;
           iLBCdec_inst^.no_of_bytes := NO_OF_BYTES_20MS;
           iLBCdec_inst^.no_of_words := NO_OF_WORDS_20MS;
           iLBCdec_inst^.state_short_len:=STATE_SHORT_LEN_20MS;
           { ULP init }
           iLBCdec_inst^.ULP_inst:=@ULP_20msTbl;
       end
       else 
       begin
           result:=2;
           exit;//exit(2);
       end;

       fillchar(iLBCdec_inst^.syntMem,LPC_FILTERORDER*sizeof(real), 0);
       move( lsfmeanTbl[0],iLBCdec_inst^.lsfdeqold,LPC_FILTERORDER*sizeof(real));

       fillchar(iLBCdec_inst^.old_syntdenum,((LPC_FILTERORDER + 1)*NSUB_MAX)*sizeof(real), 0);
       for i:=0 to NSUB_MAX-1 do
           iLBCdec_inst^.old_syntdenum[i*(LPC_FILTERORDER+1)]:=1.0;

       iLBCdec_inst^.last_lag := 20;

       iLBCdec_inst^.prevLag := 120;
       iLBCdec_inst^.per := 0.0;
       iLBCdec_inst^.consPLICount := 0;
       iLBCdec_inst^.prevPLI := 0;
       iLBCdec_inst^.prevLpc[0] := 1.0;
       fillchar(iLBCdec_inst^.prevLpc[1],LPC_FILTERORDER*sizeof(real),0);
       fillchar(iLBCdec_inst^.prevResidual, BLOCKL_MAX*sizeof(real), 0);
       iLBCdec_inst^.seed:=777;

       fillchar(iLBCdec_inst^.hpomem, 4*sizeof(real), 0);

       iLBCdec_inst^.use_enhancer := use_enhancer;
       fillchar(iLBCdec_inst^.enh_buf, ENH_BUFL*sizeof(real), 0);
       for i:=0 to ENH_NBLOCKS_TOT-1 do
           iLBCdec_inst^.enh_period[i]:=40.0;

       iLBCdec_inst^.prev_enh_pl := 0;

       result:= Smallint(iLBCdec_inst^.blockl);
   end;

   {----------------------------------------------------------------*
    *  frame residual decoder function (subrutine to iLBC_decode)
    *---------------------------------------------------------------}

   procedure Decode(
       iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) the decoder state
                                                structure }
       decresidual:pareal;             { (o) decoded residual frame }
       start:integer;                      { (i) location of start
                                              state }
       idxForMax:integer;                  { (i) codebook index for the
                                              maximum value }
       idxVec:painteger;                { (i) codebook indexes for the
                                              samples  in the start
                                              state }
       syntdenum:pareal;               { (i) the decoded synthesis
                                              filter coefficients }
       cb_index:painteger;                  { (i) the indexes for the
                                              adaptive codebook }
       gain_index:painteger;            { (i) the indexes for the
                                              corresponding gains }
       extra_cb_index:painteger;        { (i) the indexes for the
                                              adaptive codebook part
                                              of start state }
       extra_gain_index:painteger;          { (i) the indexes for the
                                              corresponding gains }
       state_first:integer                 { (i) 1 if non adaptive part
                                              of start state comes
                                              first 0 if that part
                                              comes last }
   );
   var
       reverseDecresidual:array [0..BLOCKL_MAX-1] of real;
       mem:array [0..CB_MEML-1] of real;
       k, meml_gotten, Nfor, Nback, i:integer;
       diff, start_pos:integer;
       subcount, subframe:integer;
   begin

       diff := STATE_LEN - iLBCdec_inst^.state_short_len;

       if (state_first = 1) then
       begin
           start_pos := (start-1)*SUBL;
       end
       else 
       begin
           start_pos := (start-1)*SUBL + diff;
       end;

       { decode scalar part of start state }

       StateConstructW(idxForMax, idxVec,@syntdenum[(start-1)*(LPC_FILTERORDER+1)],@decresidual[start_pos], iLBCdec_inst^.state_short_len);


       if (state_first<>0) then
       begin { put adaptive part in the end }

           { setup memory }

           fillchar(mem,(CB_MEML-iLBCdec_inst^.state_short_len)*sizeof(real), 0);
           move(decresidual[start_pos],mem[CB_MEML-iLBCdec_inst^.state_short_len],iLBCdec_inst^.state_short_len*sizeof(real));

           { construct decoded vector }

           iCBConstruct(@decresidual[start_pos+iLBCdec_inst^.state_short_len],
               extra_cb_index, extra_gain_index, @mem[CB_MEML-stMemLTbl],
               stMemLTbl, diff, CB_NSTAGES);

       end
       else 
       begin{ put adaptive part in the beginning }

           { create reversed vectors for prediction }

           for k:=0 to diff-1 do
           begin
               reverseDecresidual[k] :=decresidual[(start+1)*SUBL-1-(k+iLBCdec_inst^.state_short_len)];
           end;

           { setup memory }
           meml_gotten := iLBCdec_inst^.state_short_len;
           k:=0;
           while k<meml_gotten do
           begin
               mem[CB_MEML-1-k] := decresidual[start_pos + k];
               inc(k);
           end;
           fillchar(mem, (CB_MEML-k)*sizeof(real), 0);

           { construct decoded vector }

           iCBConstruct(@reverseDecresidual[0], extra_cb_index,
               extra_gain_index, @mem[CB_MEML-stMemLTbl], stMemLTbl,
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

       Nfor := iLBCdec_inst^.nsub-start-1;

       if ( Nfor > 0 ) then
       begin
           { setup memory }
           fillchar(mem, (CB_MEML-STATE_LEN)*sizeof(real), 0);
           move( decresidual[(start-1)*SUBL],mem[CB_MEML-STATE_LEN],STATE_LEN*sizeof(real));

           { loop over sub-frames to encode }

           for subframe:=0 to Nfor-1 do
           begin
               { construct decoded vector }

               iCBConstruct(@decresidual[(start+1+subframe)*SUBL], @cb_index[subcount*CB_NSTAGES],@gain_index[subcount*CB_NSTAGES],@mem[CB_MEML-memLfTbl[subcount]], memLfTbl[subcount], SUBL, CB_NSTAGES);

               { update memory }

               move( mem[SUBL],mem[0], (CB_MEML-SUBL)*sizeof(real));
               move(decresidual[(start+1+subframe)*SUBL],mem[CB_MEML-SUBL],SUBL*sizeof(real));
               inc(subcount);
           end;

       end;

       { backward prediction of sub-frames }

       Nback := start-1;

       if ( Nback > 0 ) then
       begin

           { setup memory }
           meml_gotten := SUBL*(iLBCdec_inst^.nsub+1-start);

           if ( meml_gotten > CB_MEML ) then
           begin
               meml_gotten:=CB_MEML;
           end;
           k:=0;
           while k<meml_gotten do
           begin
               mem[CB_MEML-1-k] := decresidual[(start-1)*SUBL + k];
               inc(k);
           end;
           fillchar(mem, (CB_MEML-k)*sizeof(real), 0);

           { loop over subframes to decode }

           for subframe:=0 to Nback-1 do
           begin
               { construct decoded vector }

               iCBConstruct(@reverseDecresidual[subframe*SUBL],
                   @cb_index[subcount*CB_NSTAGES],
                   @gain_index[subcount*CB_NSTAGES],
                   @mem[CB_MEML-memLfTbl[subcount]], memLfTbl[subcount],
                   SUBL, CB_NSTAGES);

               { update memory }

               move( mem[SUBL],mem[0], (CB_MEML-SUBL)*sizeof(real));
               move(reverseDecresidual[subframe*SUBL],mem[CB_MEML-SUBL],SUBL*sizeof(real));
               inc(subcount);
           end;

           { get decoded residual from reversed vector }

           for i:=0 to SUBL*Nback-1 do
               decresidual[SUBL*Nback - i - 1] := reverseDecresidual[i];
       end;
   end;

   {----------------------------------------------------------------*
    *  main decoder function
    *---------------------------------------------------------------}

   procedure iLBC_decode(
       decblock:pareal;            { (o) decoded signal block }
       bytes:pchar;           { (i) encoded signal bits }
       iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) the decoder state
                                                structure }
       mode:integer                    { (i) 0: bad packet, PLC,
                                              1: normal }
   );
   var
       data:array [0..BLOCKL_MAX-1] of real;
       lsfdeq:array [0..LPC_FILTERORDER*LPC_N_MAX-1] of real;
       PLCresidual:array [0..BLOCKL_MAX-1] of real;
       PLClpc:array [0..LPC_FILTERORDER] of real;
       zeros:array [0..BLOCKL_MAX-1] of real;
       one:array [0..LPC_FILTERORDER] of real;
       k, i, start, idxForMax, pos, lastpart, ulp:integer;
       lag, ilag:integer;
       cc, maxcc:real;
       idxVec:array [0..STATE_LEN-1] of integer;
       check:integer;
       gain_index:array [0..NASUB_MAX*CB_NSTAGES-1] of integer;
       extra_gain_index:array [0..CB_NSTAGES-1] of integer;
       cb_index:array [0..CB_NSTAGES*NASUB_MAX-1] of integer;
       extra_cb_index:array [0..CB_NSTAGES-1] of integer;
       lsf_i:array [0..LSF_NSPLIT*LPC_N_MAX-1] of integer;
       state_first:integer;
       last_bit:integer;
       pbytes:pchar;
       weightdenum:array [0..(LPC_FILTERORDER + 1)*NSUB_MAX-1] of real;
       order_plus_one:integer;
       syntdenum:array [0..NSUB_MAX*(LPC_FILTERORDER+1)-1] of real;
       decresidual:array [0..BLOCKL_MAX-1] of real;
   begin
       if (mode>0) then
       begin { the data are good }

           { decode data }

           pbytes:=bytes;
           pos:=0;

           { Set everything to zero before decoding }

           for k:=0 to LSF_NSPLIT*LPC_N_MAX-1 do
           begin
               lsf_i[k]:=0;
           end;
           start:=0;
           state_first:=0;
           idxForMax:=0;
           for k:=0 to iLBCdec_inst^.state_short_len-1 do
           begin
               idxVec[k]:=0;
           end;
           for k:=0 to CB_NSTAGES-1 do
           begin
               extra_cb_index[k]:=0;
               extra_gain_index[k]:=0;
           end;
           for i:=0 to iLBCdec_inst^.nasub-1 do
           begin
               for k:=0 to CB_NSTAGES-1 do
               begin
                   cb_index[i*CB_NSTAGES+k]:=0;
                   gain_index[i*CB_NSTAGES+k]:=0;
               end;
           end;

           { loop over ULP classes }

           for ulp:=0 to 2 do
           begin
               { LSF }
               for k:=0 to LSF_NSPLIT*iLBCdec_inst^.lpc_n-1 do
               begin
                   unpack( @pbytes, @lastpart,iLBCdec_inst^.ULP_inst^.lsf_bits[k][ulp], @pos);
                   packcombine(@lsf_i[k], lastpart,iLBCdec_inst^.ULP_inst^.lsf_bits[k][ulp]);
               end;

               { Start block info }

               unpack( @pbytes, @lastpart,iLBCdec_inst^.ULP_inst^.start_bits[ulp], @pos);
               packcombine(@start, lastpart,iLBCdec_inst^.ULP_inst^.start_bits[ulp]);

               unpack( @pbytes, @lastpart,iLBCdec_inst^.ULP_inst^.startfirst_bits[ulp], @pos);
               packcombine(@state_first, lastpart,iLBCdec_inst^.ULP_inst^.startfirst_bits[ulp]);

               unpack( @pbytes, @lastpart,iLBCdec_inst^.ULP_inst^.scale_bits[ulp], @pos);
               packcombine(@idxForMax, lastpart,iLBCdec_inst^.ULP_inst^.scale_bits[ulp]);

               for k:=0 to iLBCdec_inst^.state_short_len-1 do
               begin
                   unpack( @pbytes, @lastpart,iLBCdec_inst^.ULP_inst^.state_bits[ulp], @pos);
                   packcombine(@idxVec[k], lastpart,iLBCdec_inst^.ULP_inst^.state_bits[ulp]);
               end;

               { 23/22 (20ms/30ms) sample block }

               for k:=0 to CB_NSTAGES-1 do
               begin
                   unpack( @pbytes, @lastpart,iLBCdec_inst^.ULP_inst^.extra_cb_index[k][ulp],@pos);
                   packcombine(@extra_cb_index[k], lastpart,iLBCdec_inst^.ULP_inst^.extra_cb_index[k][ulp]);
               end;
               for k:=0 to CB_NSTAGES-1 do
               begin
                   unpack( @pbytes, @lastpart,iLBCdec_inst^.ULP_inst^.extra_cb_gain[k][ulp],@pos);
                   packcombine(@extra_gain_index[k], lastpart,iLBCdec_inst^.ULP_inst^.extra_cb_gain[k][ulp]);
               end;

               { The two/four (20ms/30ms) 40 sample sub-blocks }

               for i:=0 to iLBCdec_inst^.nasub-1 do
               begin
                   for k:=0 to CB_NSTAGES-1 do
                   begin
                       unpack( @pbytes, @lastpart,iLBCdec_inst^.ULP_inst^.cb_index[i][k][ulp],@pos);
                       packcombine(@cb_index[i*CB_NSTAGES+k], lastpart,iLBCdec_inst^.ULP_inst^.cb_index[i][k][ulp]);
                   end;
               end;

               for i:=0 to iLBCdec_inst^.nasub-1 do
               begin
                   for k:=0 to CB_NSTAGES-1 do
                   begin
                       unpack( @pbytes, @lastpart,iLBCdec_inst^.ULP_inst^.cb_gain[i][k][ulp],@pos);
                       packcombine(@gain_index[i*CB_NSTAGES+k], lastpart,iLBCdec_inst^.ULP_inst^.cb_gain[i][k][ulp]);
                   end;
               end;
           end;
           { Extract last bit. If it is 1 this indicates an
              empty/lost frame }
           unpack( @pbytes, @last_bit, 1, @pos);

           { Check for bit errors or empty/lost frames }
           if (start<1) then
               mode := 0;
           if (iLBCdec_inst^.mode=20)  and  (start>3) then
               mode := 0;
           if (iLBCdec_inst^.mode=30)  and  (start>5) then
               mode := 0;
           if (last_bit=1) then
               mode := 0;

           if (mode=1) then
           begin { No bit errors was detected,
                             continue decoding }
               { adjust index }
               index_conv_dec(@cb_index);
               { decode the lsf }

               SimplelsfDEQ(@lsfdeq, @lsf_i, iLBCdec_inst^.lpc_n);
               check:=LSF_check(@lsfdeq, LPC_FILTERORDER,
                   iLBCdec_inst^.lpc_n);
               DecoderInterpolateLSF(@syntdenum, @weightdenum,
                   @lsfdeq, LPC_FILTERORDER, iLBCdec_inst);

               Decode(iLBCdec_inst, @decresidual, start, idxForMax,
                   @idxVec, @syntdenum, @cb_index, @gain_index,
                   @extra_cb_index, @extra_gain_index,
                   state_first);

               { preparing the plc for a future loss! }

               doThePLC(@PLCresidual, @PLClpc, 0, @decresidual,
                   @syntdenum [(LPC_FILTERORDER + 1)*(iLBCdec_inst^.nsub - 1)],
                   iLBCdec_inst^.last_lag, iLBCdec_inst);

               move( PLCresidual[0],decresidual[0],iLBCdec_inst^.blockl*sizeof(real));
           end;

       end;

       if (mode = 0) then
       begin
           { the data is bad (either a PLC call
            * was made or a severe bit error was detected)
            }

           { packet loss conceal }

           fillchar(zeros, BLOCKL_MAX*sizeof(real), 0);

           one[0] := 1;
           fillchar(one[1], LPC_FILTERORDER*sizeof(real), 0);

           start:=0;

           doThePLC(@PLCresidual, @PLClpc, 1, @zeros, @one,iLBCdec_inst^.last_lag, iLBCdec_inst);
           move( PLCresidual[0],decresidual[0],iLBCdec_inst^.blockl*sizeof(real));

           order_plus_one := LPC_FILTERORDER + 1;
           for i := 0 to  iLBCdec_inst^.nsub-1 do
           begin
               move( PLClpc[0],syntdenum[(i*order_plus_one)],order_plus_one*sizeof(real));
           end;
       end;

       if (iLBCdec_inst^.use_enhancer = 1) then
       begin
           { post filtering }

           iLBCdec_inst^.last_lag :=enhancerInterface(@data, @decresidual, iLBCdec_inst);

           { synthesis filtering }

           if (iLBCdec_inst^.mode=20) then
           begin
               { Enhancer has 40 samples delay }
               i:=0;
               syntFilter(@data [ i*SUBL],@iLBCdec_inst^.old_syntdenum[(i+iLBCdec_inst^.nsub-1)*(LPC_FILTERORDER+1)],SUBL, @iLBCdec_inst^.syntMem);
               for i:=1 to iLBCdec_inst^.nsub-1 do
               begin
                   syntFilter(@data [ i*SUBL],@syntdenum [ (i-1)*(LPC_FILTERORDER+1)],SUBL, @iLBCdec_inst^.syntMem);
               end;
           end
           else 
           if (iLBCdec_inst^.mode=30) then
           begin
               { Enhancer has 80 samples delay }
               for i:=0 to 1 do
               begin
                   syntFilter(@data [ i*SUBL],@iLBCdec_inst^.old_syntdenum [(i+iLBCdec_inst^.nsub-2)*(LPC_FILTERORDER+1)],SUBL, @iLBCdec_inst^.syntMem);
               end;
               for i:=2 to iLBCdec_inst^.nsub-1 do
               begin
                   syntFilter(@data [ i*SUBL],@syntdenum [ (i-2)*(LPC_FILTERORDER+1)], SUBL,@iLBCdec_inst^.syntMem);
               end;
           end;

       end
       else
       begin
           { Find last lag }
           lag := 20;
           maxcc := xCorrCoef(@decresidual[BLOCKL_MAX-ENH_BLOCKL],@decresidual[BLOCKL_MAX-ENH_BLOCKL-lag], ENH_BLOCKL);

           for ilag:=21 to 119 do
           begin
               cc := xCorrCoef(@decresidual[BLOCKL_MAX-ENH_BLOCKL],@decresidual[BLOCKL_MAX-ENH_BLOCKL-ilag],ENH_BLOCKL);

               if (cc > maxcc) then
               begin
                   maxcc := cc;
                   lag := ilag;
               end;
           end;
           iLBCdec_inst^.last_lag := lag;

           { copy data and run synthesis filter }

           move( decresidual[0],data[0],iLBCdec_inst^.blockl*sizeof(real));
           for i:=0 to iLBCdec_inst^.nsub-1 do
           begin
               syntFilter(@data[i*SUBL],@syntdenum [ i*(LPC_FILTERORDER+1)], SUBL,@iLBCdec_inst^.syntMem);
           end;
       end;

       { high pass filtering on output if desired, otherwise
          copy to out }

       hpOutput(@data, iLBCdec_inst^.blockl,decblock,@iLBCdec_inst^.hpomem);

       { memcpy(decblock,data,iLBCdec_inst^.blockl*sizeof(real));}

       move( syntdenum[0],iLBCdec_inst^.old_syntdenum,iLBCdec_inst^.nsub*(LPC_FILTERORDER+1)*sizeof(real));

       iLBCdec_inst^.prev_enh_pl:=0;

       if (mode=0) then
       begin { PLC was used }
           iLBCdec_inst^.prev_enh_pl:=1;
       end;
   end;
end.