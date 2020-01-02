{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit gainquants;

{$MODE Delphi}

interface
uses constants,filter,C2Delphi_header;

   {----------------------------------------------------------------*
    *  quantizer for the gain in the gain-shape coding of residual
    *---------------------------------------------------------------}
   function gainquant({ (o) quantized gain value }
       nin:real;       { (i) gain value }
       maxIn:real;{ (i) maximum of gain value }
       cblen:integer;      { (i) number of quantization indices }
       index:pinteger      { (o) quantization index }
   ):real;
   function gaindequant(  { (o) quantized gain value }
       index:integer;      { (i) quantization index }
       maxIn:real;{ (i) maximum of unquantized gain }
       cblen:integer       { (i) number of quantization indices }
   ):real;
implementation
   function gainquant({ (o) quantized gain value }
       nin:real;       { (i) gain value }
       maxIn:real;{ (i) maximum of gain value }
       cblen:integer;      { (i) number of quantization indices }
       index:pinteger      { (o) quantization index }
   ):real;
   var
       i, tindex:integer;
       minmeasure,measure, scale:real;
       cb:pareal;
   begin
       { ensure a lower bound on the scaling factor }

       scale:=maxIn;

       if (scale<0.1) then
       begin
           scale:=0.1;
       end;

       { select the quantization table }

       if (cblen = 8) then
       begin
           cb := @gain_sq3Tbl;
       end
       else 
       if (cblen = 16) then
       begin
           cb := @gain_sq4Tbl;
       end
       else
       begin
           cb := @gain_sq5Tbl;
       end;

       { select the best index in the quantization table }

       minmeasure:=10000000.0;
       tindex:=0;
       for i:=0 to cblen-1 do
       begin
           measure:=(nin-scale*cb[i])*(nin-scale*cb[i]);

           if (measure<minmeasure) then
           begin
               tindex:=i;
               minmeasure:=measure;
           end;
       end;
       index^:=tindex;

       { return the quantized value }

       result:=scale*cb[tindex];
   end;

   {----------------------------------------------------------------*
    *  decoder for quantized gains in the gain-shape coding of
    *  residual
    *---------------------------------------------------------------}

   function gaindequant(  { (o) quantized gain value }
       index:integer;      { (i) quantization index }
       maxIn:real;{ (i) maximum of unquantized gain }
       cblen:integer       { (i) number of quantization indices }
   ):real;
   var
   	scale:real;
   begin
       { obtain correct scale factor }

       scale:=abs(maxIn);

       if (scale<0.1) then
       begin
           scale:=0.1;
       end;

       { select the quantization table and return the decoded value }

       if (cblen=8) then
       begin
           result:=scale*gain_sq3Tbl[index];
           exit;
       end
       else
       if (cblen=16) then
       begin
           result:= scale*gain_sq4Tbl[index];
           exit;
       end
       else 
       if (cblen=32) then
       begin
           result:=scale*gain_sq5Tbl[index];
           exit;
       end;

       result:=0.0;
   end;
end.
