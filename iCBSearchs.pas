{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit iCBSearchs;

{$MODE Delphi}

interface
uses iLBC_define,gainquants,createCB,filter,constants,C2Delphi_header;

   {----------------------------------------------------------------*
    *  Search routine for codebook encoding and gain quantization.
    *---------------------------------------------------------------}
   procedure iCBSearch(
       iLBCenc_inst:piLBC_Enc_Inst_t;
                           { (i) the encoder state structure }
       index:painteger;         { (o) Codebook indices }
       gain_index:painteger;{ (o) Gain quantization indices }
       intarget:pareal;{ (i) Target vector for encoding }
       mem:pareal;         { (i) Buffer for codebook construction }
       lMem:integer;           { (i) Length of buffer }
       lTarget:integer;    { (i) Length of vector }
       nStages:integer;    { (i) Number of codebook stages }
       weightDenum:pareal; { (i) weighting filter coefficients }
       weightState:pareal; { (i) weighting filter state }
       block:integer           { (i) the sub-block number }
   );
implementation
   procedure iCBSearch(
       iLBCenc_inst:piLBC_Enc_Inst_t;
                           { (i) the encoder state structure }
       index:painteger;         { (o) Codebook indices }
       gain_index:painteger;{ (o) Gain quantization indices }
       intarget:pareal;{ (i) Target vector for encoding }
       mem:pareal;         { (i) Buffer for codebook construction }
       lMem:integer;           { (i) Length of buffer }
       lTarget:integer;    { (i) Length of vector }
       nStages:integer;    { (i) Number of codebook stages }
       weightDenum:pareal; { (i) weighting filter coefficients }
       weightState:pareal; { (i) weighting filter state }
       block:integer           { (i) the sub-block number }
   );
   var
       i, j, icount, stage, best_index, range, counter:integer;
       max_measure, gain, measure, crossDot, ftmp:real;
       gains:array [0..CB_NSTAGES-1] of real;
       target:array [0..SUBL-1] of real;
       base_index, sInd, eInd, base_size:integer;
       sIndAug, eIndAug:integer;
       buf:array [0..CB_MEML+SUBL+2*LPC_FILTERORDER-1] of real;
       invenergy:array [0..CB_EXPAND*128-1] of real;
       energy:array [0..CB_EXPAND*128-1] of real;
       pp, ppi,ppo,ppe,ppe1:^real;
       cbvectors:array [0..CB_MEML-1] of real;
       tene, cene:real; 
       cvec:array [0..SUBL-1] of real;
       aug_vec:array [0..SUBL-1] of real;
			 filterno, position:integer;
   begin
       //sIndAug:=0;
       eIndAug:=0;
       ppi:=nil;
       ppo:=nil;
       ppe:=nil;
       
       fillchar(cvec,SUBL*sizeof(real),0);

       { Determine size of codebook sections }

       base_size:=lMem-lTarget+1;

       if (lTarget=SUBL) then
       begin
           base_size:=lMem-lTarget+1+lTarget div 2;
       end;

       { setup buffer for weighting }

       move(weightState[0],buf[0],sizeof(real)*LPC_FILTERORDER);
       move(mem[0],buf[LPC_FILTERORDER],lMem*sizeof(real));
       move(intarget[0],buf[LPC_FILTERORDER+lMem],lTarget*sizeof(real));

       { weighting }

       AllPoleFilter(@buf[LPC_FILTERORDER], weightDenum,lMem+lTarget, LPC_FILTERORDER);

       { Construct the codebook and target needed }

       move( buf[LPC_FILTERORDER+lMem],target[0], lTarget*sizeof(real));

       tene:=0.0;

       for i:=0 to lTarget-1 do
       begin
           tene:=tene+target[i]*target[i];
       end;

       { Prepare search over one more codebook section. This section
          is created by filtering the original buffer with a filter. }

       filteredCBvecs(@cbvectors, @buf[LPC_FILTERORDER], lMem);

       { The Main Loop over stages }

       for stage:=0 to nStages-1 do
       begin
           range := search_rangeTbl[block][stage];
           { initialize search measure }

           max_measure := -10000000.0;
           gain := 0.0;
           best_index := 0;

           { Compute cross dot product between the target
              and the CB memory }

           crossDot:=0.0;
           pp:=@buf[LPC_FILTERORDER+lMem-lTarget];
           for j:=0 to lTarget-1 do
           begin
               crossDot :=crossDot + target[j]*(pp^);
               inc(pp);
           end;

           if (stage=0) then
           begin
               { Calculate energy in the first block of
                 'lTarget' samples. }
               ppe := @energy[0];
               ppi := @buf[LPC_FILTERORDER+lMem-lTarget-1];
               ppo := @buf[LPC_FILTERORDER+lMem-1];

               ppe^:=0.0;
               pp:=@buf[LPC_FILTERORDER+lMem-lTarget];
               for j:=0 to lTarget-1 do
               begin
                   ppe^:=ppe^+(pp^)*(pp^);
                   inc(pp);
               end;

               if (ppe^>0.0) then
               begin
                   invenergy[0] :=  1.0 / (ppe^ + EPS);
               end
               else 
               begin
                   invenergy[0] := 0.0;
               end;
               inc(ppe);

               measure:=-10000000.0;

               if (crossDot > 0.0) then
               begin
                  measure := crossDot*crossDot*invenergy[0];
               end;
           end
           else 
           begin
               measure := crossDot*crossDot*invenergy[0];
           end;

           { check if measure is better }
           ftmp := crossDot*invenergy[0];

           if ((measure>max_measure)  and  (abs(ftmp)<CB_MAXGAIN)) then
           begin
               best_index := 0;
               max_measure := measure;
               gain := ftmp;
           end;

           { loop over the main first codebook section,
              full search }

           for icount:=1 to range-1 do
           begin

               { calculate measure }

               crossDot:=0.0;
               pp := @buf[LPC_FILTERORDER+lMem-lTarget-icount];

               for j:=0 to lTarget-1 do
               begin
                   crossDot :=crossDot + target[j]*(pp^);
                   inc(pp);
               end;

               if (stage=0) then
               begin
                   ppe^ := energy[icount-1] + (ppi^)*(ppi^) -(ppo^)*(ppo^);
                   inc(ppe);
                   dec(ppo);
                   dec(ppi);

                   if (energy[icount]>0.0) then
                   begin
                       invenergy[icount] :=1.0/(energy[icount]+EPS);
                   end
                   else 
                   begin
                       invenergy[icount] := 0.0;
                   end;

                   measure:=-10000000.0;

                   if (crossDot > 0.0) then
                   begin
                       measure := crossDot*crossDot*invenergy[icount];
                   end;
               end
               else 
               begin
                   measure := crossDot*crossDot*invenergy[icount];
               end;

               { check if measure is better }
               ftmp := crossDot*invenergy[icount];

               if ((measure>max_measure)  and  (abs(ftmp)<CB_MAXGAIN)) then
               begin
                   best_index := icount;
                   max_measure := measure;
                   gain := ftmp;
               end;
           end;

           { Loop over augmented part in the first codebook
            * section, full search.
            * The vectors are interpolated.
            }

           if (lTarget=SUBL) then
           begin

               { Search for best possible cb vector and
                  compute the CB-vectors' energy. }
               searchAugmentedCB(20, 39, stage, base_size-lTarget div 2,
                   @target, @buf[LPC_FILTERORDER+lMem],
                   @max_measure, @best_index, @gain, @energy,
                   @invenergy);
           end;

           { set search range for following codebook sections }

           base_index:=best_index;

           { unrestricted search }

           if (CB_RESRANGE = -1) then
           begin
               sInd:=0;
               eInd:=range-1;
               sIndAug:=20;
               eIndAug:=39;
           end

           { restricted search around best index from first
           codebook section }

           else 
           begin
               { Initialize search indices }
               sIndAug:=0;
               eIndAug:=0;
               sInd:=base_index-CB_RESRANGE div 2;
               eInd:=sInd+CB_RESRANGE;

               if (lTarget=SUBL) then
               begin
                   if (sInd<0) then
                   begin
                       sIndAug := 40 + sInd;
                       eIndAug := 39;
                       sInd:=0;
                   end
                   else 
                   if ( base_index < (base_size-20) ) then
                   begin
                       if (eInd > range) then
                       begin
                           sInd :=sInd - (eInd-range);
                           eInd := range;
                       end;
                   end 
                   else 
                   begin { base_index >= (base_size-20) }

                       if (sInd < (base_size-20)) then
                       begin
                           sIndAug := 20;
                           sInd := 0;
                           eInd := 0;
                           eIndAug := 19 + CB_RESRANGE;
                           if(eIndAug > 39) then
                           begin
                               eInd := eIndAug-39;
                               eIndAug := 39;
                           end;
                       end 
                       else 
                       begin
                           sIndAug := 20 + sInd - (base_size-20);
                           eIndAug := 39;
                           sInd := 0;
                           eInd := CB_RESRANGE - (eIndAug-sIndAug+1);
                       end;
                   end;

               end 
               else 
               begin { lTarget := 22 or 23 }
                   if (sInd < 0) then
                   begin
                       eInd :=eInd - sInd;
                       sInd := 0;
                   end;
                   if(eInd > range) then
                   begin
                       sInd :=sInd - (eInd - range);
                       eInd := range;
                   end;
               end;
           end;

           { search of higher codebook section }

           { index search range }
           counter := sInd;
           sInd :=sInd + base_size;
           eInd :=eInd + base_size;

           if (stage=0) then
           begin
               ppe := @energy[base_size];
               ppe^:=0.0;

               pp:=@cbvectors[lMem-lTarget];
               for j:=0 to lTarget-1 do
               begin
                   ppe^:=ppe^ + (pp^)*(pp^);
                   inc(pp);
               end;

               ppi := @cbvectors [ lMem - 1 - lTarget];
               ppo := @cbvectors [ lMem - 1];
               //auxulli
               ppe1:=ppe;
               inc(ppe1);
               for j:=0 to (range-2) do
               begin
                   ppe1^ := ppe^ + (ppi^)*(ppi^) - (ppo^)*(ppo^);
                   dec(ppo);
                   dec(ppi);
                   inc(ppe);
                   inc(ppe1);
               end;
           end;

           { loop over search range }

           for icount:=sInd to eInd-1 do
           begin

               { calculate measure }

               crossDot:=0.0;
               pp:=@cbvectors [ lMem - (counter) - lTarget];
               inc(counter);
               for j:=0 to lTarget-1 do
               begin
                   crossDot :=crossDot + target[j]*(pp^);
                   inc(pp);
               end;

               if (energy[icount]>0.0) then
               begin
                   invenergy[icount] :=1.0/(energy[icount]+EPS);
               end
               else
               begin
                   invenergy[icount] :=0.0;
               end;

               if (stage=0) then
               begin
                   measure:=-10000000.0;

                   if (crossDot > 0.0) then
                   begin
                       measure := crossDot*crossDot*invenergy[icount];
                   end;
               end
               else
               begin
                   measure := crossDot*crossDot*invenergy[icount];
               end;

               { check if measure is better }
               ftmp := crossDot*invenergy[icount];

               if ((measure>max_measure)  and  (abs(ftmp)<CB_MAXGAIN)) then
               begin
                   best_index := icount;
                   max_measure := measure;
                   gain := ftmp;
               end;
           end;

           { Search the augmented CB inside the limited range. }

           if ((lTarget=SUBL) and (sIndAug<>0)) then
           begin
               searchAugmentedCB(sIndAug, eIndAug, stage,
                   2*base_size-20, @target, @cbvectors[lMem],
                   @max_measure, @best_index, @gain, @energy,
                   @invenergy);
           end;

           { record best index }

           index[stage] := best_index;

           { gain quantization }

           if (stage=0) then
           begin
               if (gain<0.0) then
               begin
                   gain := 0.0;
               end;

               if (gain>CB_MAXGAIN) then
               begin
                   gain := CB_MAXGAIN;
               end;
               gain := gainquant(gain, 1.0, 32, @gain_index[stage]);
           end
           else 
           begin
               if (stage=1) then
               begin
                   gain := gainquant(gain, abs(gains[stage-1]),
                       16, @gain_index[stage]);
               end 
               else 
               begin
                   gain := gainquant(gain, abs(gains[stage-1]),
                       8, @gain_index[stage]);
               end;
           end;

           { Extract the best (according to measure)
              codebook vector }

           if (lTarget=(STATE_LEN-iLBCenc_inst^.state_short_len)) then
           begin
               if (index[stage]<base_size) then
               begin
                   pp:=@buf[LPC_FILTERORDER+lMem-lTarget-index[stage]];
               end 
               else 
               begin
                   pp:=@cbvectors[lMem-lTarget-index[stage]+base_size];
               end;
           end 
           else 
           begin
               if (index[stage]<base_size) then
               begin
                   if (index[stage]<(base_size-20)) then
                   begin
                       pp:=@buf[LPC_FILTERORDER+lMem-lTarget-index[stage]];
                   end
                   else 
                   begin
                       createAugmentedVec(index[stage]-base_size+40,@buf[LPC_FILTERORDER+lMem],@aug_vec);
                       pp:=@aug_vec;
                   end;
               end 
               else 
               begin
                   filterno:=index[stage] div base_size;
                   position:=index[stage]-filterno*base_size;

                   if (position<(base_size-20)) then
                   begin
                       pp:=@cbvectors[filterno*lMem-lTarget-index[stage]+filterno*base_size];
                   end 
                   else 
                   begin
                       createAugmentedVec(index[stage]-(filterno+1)*base_size+40,@cbvectors[filterno*lMem],@aug_vec);
                       pp:=@aug_vec;
                   end;
               end;
           end;

           { Subtract the best codebook vector, according
              to measure, from the target vector }

           for j:=0 to lTarget-1 do
           begin
               cvec[j] :=cvec[j] + gain*(pp^);
               target[j] :=target[j] - gain*(pp^);
               inc(pp);
           end;

           { record quantized gain }

           gains[stage]:=gain;

       end;{ end of Main Loop. for (stage:=0;... }

       { Gain adjustment for energy matching }
       cene:=0.0;
       for i:=0 to lTarget-1 do
       begin
           cene:=cene+cvec[i]*cvec[i];
       end;
       j:=gain_index[0];

       for i:=gain_index[0] to 31 do
       begin
           ftmp:=cene*gain_sq5Tbl[i]*gain_sq5Tbl[i];
           if ((ftmp<(tene*gains[0]*gains[0]))  and  (gain_sq5Tbl[j]<(2.0*gains[0]))) then
           begin
               j:=i;
           end;
       end;
       gain_index[0]:=j;
   end;
end.





