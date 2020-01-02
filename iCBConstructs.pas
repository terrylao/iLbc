unit iCBConstructs;

{$MODE Delphi}

interface
uses iLBC_define,gainquants,getCBvecs,C2Delphi_header;

   {----------------------------------------------------------------*
    *  Convert the codebook indexes to make the search easier
    *---------------------------------------------------------------}
   procedure index_conv_enc(
       index:painteger          { (i/o) Codebook indexes }
   );
   procedure index_conv_dec(
       index:painteger          { (i/o) Codebook indexes }
   );
   procedure iCBConstruct(
       decvector:pareal;   { (o) Decoded vector }
       index:painteger;         { (i) Codebook indices }
       gain_index:painteger;{ (i) Gain quantization indices }
       mem:pareal;         { (i) Buffer for codevector construction }
       lMem:integer;           { (i) Length of buffer }
       veclen:integer;         { (i) Length of vector }
       nStages:integer         { (i) Number of codebook stages }
   );
implementation
   procedure index_conv_enc(
       index:painteger          { (i/o) Codebook indexes }
   );
   var
   	k:integer;
   begin
       for k:=1 to CB_NSTAGES-1 do
       begin
           if ((index[k]>=108) and (index[k]<172)) then
           begin
               index[k]:=index[k]-64
           end 
           else 
           if (index[k]>=236) then
           begin
               index[k]:=index[k]-128;
           end
           else 
           begin
               { ERROR }
           end;
       end;
   end;

   procedure index_conv_dec(
       index:painteger          { (i/o) Codebook indexes }
   );
   var
   	k:integer;
   begin
       for k:=1 to CB_NSTAGES-1 do
       begin
           if ((index[k]>=44) and (index[k]<108)) then
           begin
               index[k]:=index[k]+64;
           end 
           else 
           if ((index[k]>=108) and (index[k]<128)) then
           begin
               index[k]:=index[k]+128;
           end 
           else 
           begin
               { ERROR }
           end;
       end;
   end;

   {----------------------------------------------------------------*
    *  Construct decoded vector from codebook and gains.
    *---------------------------------------------------------------}

   procedure iCBConstruct(
       decvector:pareal;   { (o) Decoded vector }
       index:painteger;         { (i) Codebook indices }
       gain_index:painteger;{ (i) Gain quantization indices }
       mem:pareal;         { (i) Buffer for codevector construction }
       lMem:integer;           { (i) Length of buffer }
       veclen:integer;         { (i) Length of vector }
       nStages:integer         { (i) Number of codebook stages }
   );
   var
       j,k:integer;
       gain:array [0..CB_NSTAGES-1] of real;
       cbvec:array [0..SUBL] of real;
   begin

       { gain de-quantization }

       gain[0] := gaindequant(gain_index[0], 1.0, 32);
       if (nStages > 1) then
       begin
           gain[1] := gaindequant(gain_index[1],abs(gain[0]), 16);
       end;
       if (nStages > 2) then
       begin
           gain[2] := gaindequant(gain_index[2],abs(gain[1]), 8);
       end;

       { codebook vector construction and construction of
       total vector }

       getCBvec(@cbvec, mem, index[0], lMem, veclen);
       for j:=0 to veclen-1 do
       begin
           decvector[j] := gain[0]*cbvec[j];
       end;
       if (nStages > 1) then
       begin
           for k:=1 to nStages-1 do
           begin
               getCBvec(@cbvec, mem, index[k], lMem, veclen);
               for j:=0 to veclen-1 do
               begin
                   decvector[j] :=decvector[j] + gain[k]*cbvec[j];
               end;
           end;
       end;
   end;
end.