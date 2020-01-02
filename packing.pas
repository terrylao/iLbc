{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit packing;

{$MODE Delphi}

interface
uses
  iLBC_define,
  constants,
  helpfun,C2Delphi_header;

   {----------------------------------------------------------------*
    *  splitting an integer into first most significant bits and
    *  remaining least significant bits
    *---------------------------------------------------------------}
   procedure packsplit(
       index:pinteger;                 { (i) the value to split }
       firstpart:pinteger;             { (o) the value specified by most
                                          significant bits }
       rest:pinteger;                  { (o) the value specified by least
                                          significant bits }
       bitno_firstpart:integer;    { (i) number of bits in most
                                          significant part }
       bitno_total:integer             { (i) number of bits in full range
                                          of value }
   );
   procedure packcombine(
       index:pinteger;                 { (i/o) the msb value in the
                                          combined value out }
       rest:integer;                   { (i) the lsb value }
       bitno_rest:integer              { (i) the number of bits in the
                                          lsb part }
   );
   procedure dopack(
       bitstream :ppchar;  { (i/o) on entrance pointer to
                                          place in bitstream to pack
                                          new data, on exit pointer
                                          to place in bitstream to
                                          pack future data }
       index:integer;                  { (i) the value to pack }
       bitno:integer;                  { (i) the number of bits that the
                                          value will fit within }
       pos:pinteger                { (i/o) write position in the
                                          current byte }
   );
   procedure unpack(
       bitstream :ppchar;  { (i/o) on entrance pointer to
                                          place in bitstream to
                                          unpack new data from, on
                                          exit pointer to place in
                                          bitstream to unpack future
                                          data from }
       index:pinteger;                 { (o) resulting value }
       bitno:integer;                  { (i) number of bits used to
                                          represent the value }
       pos:pinteger                { (i/o) read position in the
                                          current byte }
   );
implementation
   procedure packsplit(
       index:pinteger;                 { (i) the value to split }
       firstpart:pinteger;             { (o) the value specified by most
                                          significant bits }
       rest:pinteger;                  { (o) the value specified by least
                                          significant bits }
       bitno_firstpart:integer;    { (i) number of bits in most
                                          significant part }
       bitno_total:integer             { (i) number of bits in full range
                                          of value }
   );
   var
   	bitno_rest:integer;
   begin
       bitno_rest := bitno_total-bitno_firstpart;

       firstpart^ := index^ shr (bitno_rest);
       rest^ := index^-(firstpart^ shl (bitno_rest));
   end;

   {----------------------------------------------------------------*
    *  combining a value corresponding to msb's with a value
    *  corresponding to lsb's
    *---------------------------------------------------------------}

   procedure packcombine(
       index:pinteger;                 { (i/o) the msb value in the
                                          combined value out }
       rest:integer;                   { (i) the lsb value }
       bitno_rest:integer              { (i) the number of bits in the
                                          lsb part }
   );
   begin
       index^ := index^ shl bitno_rest;
       index^ :=index^ + rest;
   end;

   {----------------------------------------------------------------*
    *  packing of bits into bitstream, i.e., vector of bytes
    *---------------------------------------------------------------}

   procedure dopack(
       bitstream :ppchar;  { (i/o) on entrance pointer to
                                          place in bitstream to pack
                                          new data, on exit pointer
                                          to place in bitstream to
                                          pack future data }
       index:integer;                  { (i) the value to pack }
       bitno:integer;                  { (i) the number of bits that the
                                          value will fit within }
       pos:pinteger                { (i/o) write position in the
                                          current byte }
   );
   var
     posLeft:integer;
   begin
       { Clear the bits before starting in a new byte }

       if (pos^=0) then
       begin
           bitstream^^:=#0;
       end;

       while (bitno>0) do
       begin

           { Jump to the next byte if end of this byte is reached}

           if (pos^=8) then
           begin
               pos^:=0;
               inc(bitstream^);//(*bitstream)++;
               bitstream^^:=#0;
           end;

           posLeft:=8-(pos^);

           { Insert index into the bitstream }

           if (bitno <= posLeft) then
           begin
               bitstream^^ := char(byte(bitstream^^) or (index shl (posLeft-bitno)));
               pos^:=pos^+bitno;
               bitno:=0;
           end
           else
           begin
               bitstream^^:=char(byte(bitstream^^) or (index shr (bitno-posLeft)));

               pos^:=8;
               index:=index-((index shr (bitno-posLeft)) shl (bitno-posLeft));

               bitno:=bitno-posLeft;
           end;
       end;
   end;

   {----------------------------------------------------------------*
    *  unpacking of bits from bitstream, i.e., vector of bytes
    *---------------------------------------------------------------}

   procedure unpack(
       bitstream :ppchar;  { (i/o) on entrance pointer to
                                          place in bitstream to
                                          unpack new data from, on
                                          exit pointer to place in
                                          bitstream to unpack future
                                          data from }
       index:pinteger;                 { (o) resulting value }
       bitno:integer;                  { (i) number of bits used to
                                          represent the value }
       pos:pinteger                { (i/o) read position in the
                                          current byte }
   );
   var
   	BitsLeft:integer;
   begin

       index^:=0;

       while (bitno>0) do
       begin
           { move forward in bitstream when the end of the
              byte is reached }

           if (pos^=8) then
           begin
               pos^:=0;
               inc(bitstream^);//(*bitstream)++;
           end;

           BitsLeft:=8-(pos^);

           { Extract bits to index }

           if (BitsLeft>=bitno) then
           begin
               index^:=index^+ (((byte(bitstream^^) shl (pos^)) and $FF) shr (8-bitno));

               pos^:=pos^+bitno;
               bitno:=0;
           end
           else
           begin

               if ((8-bitno)>0) then
               begin
                   index^:=index^+(((byte(bitstream^^) shl (pos^)) and $FF) shr (8-bitno));
                   pos^:=8;
               end
               else
               begin
                   index^:=index^+((((byte(bitstream^^) shl (pos^)) and $FF)) shl (bitno-8));
                   pos^:=8;
               end;
               bitno:=bitno-BitsLeft;
           end;
       end;
   end;
end.
