unit anaFilters;

{$MODE Delphi}

interface
uses
   iLBC_define,C2Delphi_header;

   {----------------------------------------------------------------*
    *  LP analysis filter.
    *---------------------------------------------------------------}
    procedure anaFilter(
       nIn:pareal;  { (i) Signal to be filtered }
       a:pareal;   { (i) LP parameters }
       len:integer;{ (i) Length of signal }
       nOut:pareal; { (o) Filtered signal }
       mem:pareal  { (i/o) Filter state }
   );

implementation
   procedure anaFilter(
       nIn:pareal;  { (i) Signal to be filtered }
       a:pareal;   { (i) LP parameters }
       len:integer;{ (i) Length of signal }
       nOut:pareal; { (o) Filtered signal }
       mem:pareal  { (i/o) Filter state }
   );
   var
   	i,j:integer;
    po, pi, pm, pa:^real;
   begin
       po := @nOut[0];

       { Filter first part using memory from past }
			 
			 for i:=0 to LPC_FILTERORDER-1 do
			 begin
			 	pi:=@(nIn[i]);
			 	pm:=@(mem[LPC_FILTERORDER-1]);
			 	pa:=@a[0];
			 	po^:=0.0;
			 	for j:=0 to i do
			 	begin
			 		po^:=po^+pa^*pi^;
			 		inc(pa);
			 		dec(pi);
			 	end;
			 	for j:=i+1 to LPC_FILTERORDER do
			 	begin
			 		po^:=po^+pa^*pm^;
			 		inc(pa);
			 		dec(pm)
			 	end;
			 	inc(po);
			 end;

       { Filter last part where the state is entirely
          in the input vector }
			 for i:=LPC_FILTERORDER to len-1 do
			 begin
			 	pi:=@(nIn[i]);
			 	pa:=@a[0];
			 	po^:=0.0;
			 	for j:=0 to LPC_FILTERORDER do
			 	begin
			 		po^:=po^+pa^*pi^;
			 		inc(pa);
			 		dec(pi);
			 	end;
			 	inc(po);
			 end;

       { Update state vector }
			 Move( nIn[len-LPC_FILTERORDER],mem[0], LPC_FILTERORDER*sizeof(real));   
   end;
end.

