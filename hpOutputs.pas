unit hpOutputs;

{$MODE Delphi}

interface
uses constants,C2Delphi_header;

   {----------------------------------------------------------------*
    *  Output high-pass filter
    *---------------------------------------------------------------}
   procedure hpOutput(
       nIn:pareal;  { (i) vector to filter }
       len:integer;{ (i) length of vector to filter }
       nOut:pareal; { (o) the resulting filtered vector }
       mem:pareal  { (i/o) the filter state }
   );
implementation
   procedure hpOutput(
       nIn:pareal;  { (i) vector to filter }
       len:integer;{ (i) length of vector to filter }
       nOut:pareal; { (o) the resulting filtered vector }
       mem:pareal  { (i/o) the filter state }
   );
   var
       i:integer;
       pi, po:^real;
   begin

       { all-zero section}

       pi := @nIn[0];
       po := @nOut[0];
       for i:=0 to len-1 do
       begin
           po^ := hpo_zero_coefsTbl[0] * (pi^);
           po^ :=po^ + hpo_zero_coefsTbl[1] * mem[0];
           po^ :=po^ + hpo_zero_coefsTbl[2] * mem[1];

           mem[1] := mem[0];
           mem[0] := pi^;
           inc(po);
           inc(pi);

       end;

       { all-pole section}

       po := @nOut[0];
       for i:=0 to len-1 do
       begin
           po^ :=po^ - hpo_pole_coefsTbl[1] * mem[2];
           po^ :=po^ - hpo_pole_coefsTbl[2] * mem[3];

           mem[3] := mem[2];
           mem[2] := po^;
           inc(po);
       end;
   end;
end.