{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit createCB;

{$MODE Delphi}

interface
uses iLBC_define,constants,C2Delphi_header;
   procedure filteredCBvecs(
       cbvectors:pareal;   { (o) Codebook vectors for the
                                  higher section }
       mem:pareal;         { (i) Buffer to create codebook
                                  vector from }
       lMem:integer        { (i) Length of buffer }
   );
   procedure searchAugmentedCB(
       low:integer;        { (i) Start index for the search }
       high:integer;           { (i) End index for the search }
       stage:integer;          { (i) Current stage }
       startIndex:integer;     { (i) Codebook index for the first
                                  aug vector }
       target:PAreal;      { (i) Target vector for encoding }
       buffer:pareal;      { (i) Pointer to the end of the buffer for
                                  augmented codebook construction }
       max_measure:pareal; { (i/o) Currently maximum measure }
       best_index:painteger;{ (o) Currently the best index }
       gain:pareal;    { (o) Currently the best gain }
       energy:PAreal;      { (o) Energy of augmented codebook
                                  vectors }
       invenergy:PAreal{ (o) Inv energy of augmented codebook
                                  vectors }
   );
   procedure createAugmentedVec(
       index:integer;      { (i) Index for the augmented vector
                              to be created }
       buffer:pareal;  { (i) Pointer to the end of the buffer for
                              augmented codebook construction }
       cbVec:PAreal{ (o) The construced codebook vector }
   );
implementation
   {----------------------------------------------------------------*
    *  Construct an additional codebook vector by filtering the
    *  initial codebook buffer. This vector is then used to expand
    *  the codebook with an additional section.
    *---------------------------------------------------------------}

   procedure filteredCBvecs(
       cbvectors:pareal;   { (o) Codebook vectors for the
                                  higher section }
       mem:pareal;         { (i) Buffer to create codebook
                                  vector from }
       lMem:integer        { (i) Length of buffer }
   );
   var
   	j,k:integer;
    pp, pp1:^real;
    tempbuff2:array [0..CB_MEML+CB_FILTERLEN-1] of real;
    pos:^real;
   begin
			fillchar(tempbuff2,(CB_HALFFILTERLEN-1)*sizeof(real),0);
      move( mem[0],tempbuff2[CB_HALFFILTERLEN-1], lMem*sizeof(real));
      fillchar(tempbuff2[lMem+CB_HALFFILTERLEN-1],(CB_HALFFILTERLEN+1)*sizeof(real), 0);

       { Create codebook vector for higher section by filtering }

       { do filtering }
       pos:=@cbvectors[0];
       fillchar(pos^, lMem*sizeof(real), 0);
       for k:=0 to lMem-1 do
       begin
       	pp:=@tempbuff2[k];
       	pp1:=@cbfiltersTbl[CB_FILTERLEN-1];
       	for j:=0 to CB_FILTERLEN-1 do
       	begin
       		pos^:=pos^+pp^*pp1^;
       		inc(pp);
       		dec(pp1);
       	end;
       	inc(pos);
       end;
   end;

   {----------------------------------------------------------------*
    *  Search the augmented part of the codebook to find the best
    *  measure.
    *----------------------------------------------------------------}






   procedure searchAugmentedCB(
       low:integer;        { (i) Start index for the search }
       high:integer;           { (i) End index for the search }
       stage:integer;          { (i) Current stage }
       startIndex:integer;     { (i) Codebook index for the first
                                  aug vector }
       target:PAreal;      { (i) Target vector for encoding }
       buffer:pareal;      { (i) Pointer to the end of the buffer for
                                  augmented codebook construction }
       max_measure:pareal; { (i/o) Currently maximum measure }
       best_index:painteger;{ (o) Currently the best index }
       gain:pareal;    { (o) Currently the best gain }
       energy:PAreal;      { (o) Energy of augmented codebook
                                  vectors }
       invenergy:PAreal{ (o) Inv energy of augmented codebook
                                  vectors }
   );
   var
       icount, ilow, j, tmpIndex:integer;
       pp, ppo, ppi, ppe:^real;
       crossDot, alfa:real;
       weighted, measure, nrjRecursive:real;
       ftmp:real;   
   begin
       { Compute the energy for the first (low-5)
          noninterpolated samples }
       nrjRecursive := 0.0;
       pp := @buffer[ - low + 1];
       for j:=0 to low - 6 do
       begin
       	nrjRecursive := nrjRecursive+( (pp^)*(pp^) );
       	inc(pp);
       end;
       ppe := @buffer[ - low];

			 for icount:=low to high do
			 begin

           { Index of the codebook vector used for retrieving
              energy values }
           tmpIndex := startIndex+icount-20;

           ilow := icount-4;

           { Update the energy recursively to save complexity }
           nrjRecursive := nrjRecursive + (ppe^)*(ppe^);
           dec(ppe);
           energy[tmpIndex] := nrjRecursive;

           { Compute cross dot product for the first (low-5)
              samples }





           crossDot := 0.0;
           pp := @buffer[-icount];
           for j:=0 to ilow-1 do
           begin
             crossDot := crossDot+target[j]*(pp^);
             inc(pp);
           end;

           { interpolation }
           alfa := 0.2;
           j:=0;
           ppo := @buffer[j-4];
           ppi := @buffer[-icount-4];
           for j:=ilow to icount-1 do
           begin
               weighted := (1.0-alfa)*(ppo^)+alfa*(ppi^);
               inc(ppo);
               inc(ppi);
               energy[tmpIndex] := energy[tmpIndex]+weighted*weighted;
               crossDot := crossDot+target[j]*weighted;
               alfa :=alfa+ 0.2;
           end;

           { Compute energy and cross dot product for the
              remaining samples }
           pp := @buffer[ - icount];
           for j:=icount to SUBL-1 do
           begin
               energy[tmpIndex] := energy[tmpIndex]+(pp^)*(pp^);
               crossDot :=crossDot+ target[j]*(pp^);
               inc(pp);
           end;

           if (energy[tmpIndex]>0.0) then
           begin
               invenergy[tmpIndex]:=1.0/(energy[tmpIndex]+EPS);
           end 
           else
           begin
               invenergy[tmpIndex] := 0.0;
           end;

           if (stage=0) then
           begin
               measure := -10000000.0;

               if (crossDot > 0.0) then
               begin
                   measure := crossDot*crossDot*invenergy[tmpIndex];
               end;
           end
           else
           begin
               measure := crossDot*crossDot*invenergy[tmpIndex];
           end;

           { check if measure is better }
           ftmp := crossDot*invenergy[tmpIndex];

           if ((measure>max_measure[0]) and (abs(ftmp)<CB_MAXGAIN)) then
           begin
               best_index[0] := tmpIndex;
               max_measure[0] := measure;
               gain[0] := ftmp;
           end;
     end;
   end;


   {----------------------------------------------------------------*
    *  Recreate a specific codebook vector from the augmented part.
    *
    *----------------------------------------------------------------}

   procedure createAugmentedVec(
       index:integer;      { (i) Index for the augmented vector
                              to be created }
       buffer:pareal;  { (i) Pointer to the end of the buffer for
                              augmented codebook construction }
       cbVec:PAreal{ (o) The construced codebook vector }
   );
   var
     ilow, j:integer;
     pp, ppo, ppi:^real;
     alfa, alfa1, weighted:real;
   begin
       ilow := index-5;

       { copy the first noninterpolated part }

       pp := @buffer[-index];
       move(pp^,cbVec^,sizeof(real)*index);

       { interpolation }

       alfa1 := 0.2;
       alfa := 0.0;
       j:=0;
       ppo := @buffer[j-5];
       ppi := @buffer[-index-5];
       for j:=ilow to index-1 do
       begin
           weighted := (1.0-alfa)*(ppo^)+alfa*(ppi^);
           inc(ppo);
           inc(ppi);
           cbVec[j] := weighted;
           alfa :=alfa + alfa1;
       end;

       { copy the second noninterpolated part }

       pp := @buffer[ - index];
       move(pp^,cbVec[index],sizeof(real)*(SUBL-index));
   end;


end.
