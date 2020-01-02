unit filter;

{$MODE Delphi}

interface
uses iLBC_define,C2Delphi_header;

   {----------------------------------------------------------------*
    *  all-pole filter
    *---------------------------------------------------------------}
   procedure AllPoleFilter(
       InOut:PAreal;   { (i/o) on entrance InOut[-orderCoef] to
                              InOut[-1] contain the state of the
                              filter (delayed samples). InOut[0] to
                              InOut[lengthInOut-1] contain the filter
                              input, on en exit InOut[-orderCoef] to
                              InOut[-1] is unchanged and InOut[0] to
                              InOut[lengthInOut-1] contain filtered
                              samples }
       Coef:PAreal;{ (i) filter coefficients, Coef[0] is assumed
                              to be 1.0 }
       lengthInOut:integer;{ (i) number of input/output samples }
       orderCoef:integer   { (i) number of filter coefficients }
   );
   procedure AllZeroFilter(
       nIn:PAreal;      { (i) In[0] to In[lengthInOut-1] contain
                              filter input samples }
       Coef:PAreal;{ (i) filter coefficients (Coef[0] is assumed
                              to be 1.0) }
       lengthInOut:integer;{ (i) number of input/output samples }
       orderCoef:integer;  { (i) number of filter coefficients }
       nOut:PAreal      { (i/o) on entrance Out[-orderCoef] to Out[-1]
                              contain the filter state, on exit Out[0]
                              to Out[lengthInOut-1] contain filtered
                              samples }
   );
   procedure ZeroPoleFilter(
       nIn:PAreal;      { (i) In[0] to In[lengthInOut-1] contain
                              filter input samples In[-orderCoef] to
                              In[-1] contain state of all-zero
                              section }
       ZeroCoef:pareal;{ (i) filter coefficients for all-zero
                              section (ZeroCoef[0] is assumed to
                              be 1.0) }
       PoleCoef:pareal;{ (i) filter coefficients for all-pole section
                              (ZeroCoef[0] is assumed to be 1.0) }
       lengthInOut:integer;{ (i) number of input/output samples }
       orderCoef:integer;  { (i) number of filter coefficients }
       nOut:pareal      { (i/o) on entrance Out[-orderCoef] to Out[-1]
                              contain state of all-pole section. On
                              exit Out[0] to Out[lengthInOut-1]
                              contain filtered samples }
   );
   procedure DownSample (
       nIn:PAreal;     { (i) input samples }
       Coef:pareal;   { (i) filter coefficients }
       lengthIn:integer;   { (i) number of input samples }
       state:pareal;  { (i) filter state }
       nOut:PAreal     { (o) downsampled output }
   );
implementation
   procedure AllPoleFilter(
       InOut:PAreal;   { (i/o) on entrance InOut[-orderCoef] to
                              InOut[-1] contain the state of the
                              filter (delayed samples). InOut[0] to
                              InOut[lengthInOut-1] contain the filter
                              input, on en exit InOut[-orderCoef] to
                              InOut[-1] is unchanged and InOut[0] to
                              InOut[lengthInOut-1] contain filtered
                              samples }
       Coef:PAreal;{ (i) filter coefficients, Coef[0] is assumed
                              to be 1.0 }
       lengthInOut:integer;{ (i) number of input/output samples }
       orderCoef:integer   { (i) number of filter coefficients }
   );
   var
   	n,k:integer;
   begin
       for n:=0 to lengthInOut-1 do
       begin
           for k:=1 to orderCoef do
           begin
               InOut[0] :=InOut[0] - Coef[k]*InOut[-k];
           end;
           InOut:=@InOut[1];
       end;
   end;

   {----------------------------------------------------------------*
    *  all-zero filter
    *---------------------------------------------------------------}

   procedure AllZeroFilter(
       nIn:PAreal;      { (i) In[0] to In[lengthInOut-1] contain
                              filter input samples }
       Coef:PAreal;{ (i) filter coefficients (Coef[0] is assumed
                              to be 1.0) }
       lengthInOut:integer;{ (i) number of input/output samples }
       orderCoef:integer;  { (i) number of filter coefficients }
       nOut:PAreal      { (i/o) on entrance Out[-orderCoef] to Out[-1]
                              contain the filter state, on exit Out[0]
                              to Out[lengthInOut-1] contain filtered
                              samples }
   );
   var
   	n,k:integer;
   begin
       for n:=0 to lengthInOut-1 do
       begin
           nOut[n] := Coef[0]*nIn[0];
           for k:=1 to orderCoef do
           begin
               nOut[n] :=nOut[n]+ Coef[k]*nIn[-k];
           end;
           nIn:=@nIn[1];
       end;
   end;

   {----------------------------------------------------------------*
    *  pole-zero filter
    *---------------------------------------------------------------}

   procedure ZeroPoleFilter(
       nIn:PAreal;      { (i) In[0] to In[lengthInOut-1] contain
                              filter input samples In[-orderCoef] to
                              In[-1] contain state of all-zero
                              section }
       ZeroCoef:pareal;{ (i) filter coefficients for all-zero
                              section (ZeroCoef[0] is assumed to
                              be 1.0) }
       PoleCoef:pareal;{ (i) filter coefficients for all-pole section
                              (ZeroCoef[0] is assumed to be 1.0) }
       lengthInOut:integer;{ (i) number of input/output samples }
       orderCoef:integer;  { (i) number of filter coefficients }
       nOut:pareal      { (i/o) on entrance Out[-orderCoef] to Out[-1]
                              contain state of all-pole section. On
                              exit Out[0] to Out[lengthInOut-1]
                              contain filtered samples }
   );
   begin
       AllZeroFilter(nIn,ZeroCoef,lengthInOut,orderCoef,nOut);
       AllPoleFilter(nOut,PoleCoef,lengthInOut,orderCoef);
   end;

   {----------------------------------------------------------------*
    * downsample (LP filter and decimation)
    *---------------------------------------------------------------}

   procedure DownSample (
       nIn:PAreal;     { (i) input samples }
       Coef:pareal;   { (i) filter coefficients }
       lengthIn:integer;   { (i) number of input samples }
       state:pareal;  { (i) filter state }
       nOut:PAreal     { (o) downsampled output }
   );
   var
       o:real;
       Out_ptr:^real;
       Coef_ptr, In_ptr:^real;
       state_ptr:^real;
       i, j, stop:integer;
   begin
			Out_ptr := @nOut[0];
       { LP filter and decimate at the same time }
      i := DELAY_DS;
       while ( i < lengthIn) do
       begin
           Coef_ptr := @Coef[0];
           In_ptr := @nIn[i];
           state_ptr := @state[FILTERORDER_DS-2];

           o := 0.0;
           
           if (i < FILTERORDER_DS) then
            stop := i + 1 
           else
            stop := FILTERORDER_DS;

           for j := 0 to stop-1 do
           begin
               o :=o + Coef_ptr^ * In_ptr^;
               inc(Coef_ptr);
               dec(In_ptr);
           end;
           for j := i + 1 to  FILTERORDER_DS-1 do
           begin
               o :=o + Coef_ptr^ * (state_ptr^);
               inc(Coef_ptr);
               dec(state_ptr);
           end;
           Out_ptr^ := o;     
           inc(Out_ptr);
           i:=i+FACTOR_DS;
       end;

       { Get the last part (use zeros as input for the future) }
       i:=(lengthIn+FACTOR_DS);
       while ( i<(lengthIn+DELAY_DS)) do
       begin
           o:=0.0;

           if (i<lengthIn) then
           begin
               Coef_ptr := @Coef[0];
               In_ptr := @nIn[i];
               for j:=0 to FILTERORDER_DS-1 do
               begin
                  o :=o + Coef_ptr^ * Out_ptr^;
	               inc(Coef_ptr);
	               dec(Out_ptr);
               end;
           end
           else
           begin
               Coef_ptr := @Coef[i-lengthIn];
               In_ptr := @nIn[lengthIn-1];
               for j:=0 to FILTERORDER_DS-(i-lengthIn)-1 do
               begin
                  o := o+ Coef_ptr^ * In_ptr^;
	               inc(Coef_ptr);
	               dec(In_ptr);
               end;
           end;
           Out_ptr^ := o;
           inc(Out_ptr);
           i:=i+FACTOR_DS;
       end;
   end;
end.