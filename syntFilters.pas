{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit syntFilters;

{$MODE Delphi}

interface
uses iLBC_define,C2Delphi_header;

   {----------------------------------------------------------------*
    *  LP synthesis filter.
    *---------------------------------------------------------------}

   procedure syntFilter(
       nOut:pareal;     { (i/o) Signal to be filtered }
       a:pareal;       { (i) LP parameters }
       len:integer;    { (i) Length of signal }
       mem:pareal      { (i/o) Filter state }
   );
implementation
   procedure syntFilter(
       nOut:pareal;     { (i/o) Signal to be filtered }
       a:pareal;       { (i) LP parameters }
       len:integer;    { (i) Length of signal }
       mem:pareal      { (i/o) Filter state }
   );
   var
       i, j:integer;
       po, pi, pa, pm:^real;
   begin

       po:=@nOut[0];

       { Filter first part using memory from past }

       for i:=0 to LPC_FILTERORDER-1 do
       begin
           pi:=@nOut[i-1];
           pa:=@a[1];
           pm:=@mem[LPC_FILTERORDER-1];
           for j:=1 to i do
           begin
               po^:=po^-(pa^)*(pi^);
               inc(pa);
               dec(pi);
           end;
           for j:=i+1 to LPC_FILTERORDER do
           begin
               po^:=po^-(pa^)*(pm^);
               inc(pa);
               dec(pm);
           end;
           inc(po);
       end;

       { Filter last part where the state is entirely in
          the output vector }

       for i:=LPC_FILTERORDER to len-1 do
       begin
           pi:=@nOut[i-1];
           pa:=@a[1];
           for j:=1 to LPC_FILTERORDER do
           begin
               po^:=po^-(pa^)*(pi^);
               inc(pa);
               dec(pi);
           end;
           inc(po);
       end;

       { Update state vector }
       move( nOut[len-LPC_FILTERORDER],mem[0],LPC_FILTERORDER*sizeof(real));
   end;
end.
