unit helpfun;

{$MODE Delphi}

interface
uses iLBC_define,constants,C2Delphi_header;

   {----------------------------------------------------------------*
    *  calculation of auto correlation
    *---------------------------------------------------------------}
const
	 eps   = 0.039; { 50 Hz }
   eps2  =0.0195;
   maxlsf=3.14; { 4000 Hz }
   minlsf=0.01; { 0 Hz }

   procedure autocorr(
       r:pareal;       { (o) autocorrelation vector }
       x:pareal; { (i) data vector }
       N:integer;          { (i) length of data vector }
       order:integer       { largest lag for calculated
                          autocorrelations }
   );
   procedure window(
       z:pareal;       { (o) the windowed data }
       x:pareal; { (i) the original data vector }
       y:pareal; { (i) the window }
       N:integer           { (i) length of all vectors }
   );
   procedure levdurb(
       a:pareal;       { (o) lpc coefficient vector starting
                              with 1.0 }
       k:pareal;       { (o) reflection coefficients }
       r:pareal;       { (i) autocorrelation vector }
       order:integer       { (i) order of lpc filter }
   );
   procedure interpolate(
       nout:pareal;      { (o) the interpolated vector }
       in1:pareal;     { (i) the first vector for the
                              interpolation }
       in2:pareal;     { (i) the second vector for the
                              interpolation }
       coef:real;      { (i) interpolation weights }
       length:integer      { (i) length of all vectors }
   );
   procedure bwexpand(
       nout:pareal;      { (o) the bandwidth expanded lpc
                              coefficients }
       nin:pareal;      { (i) the lpc coefficients before bandwidth
                              expansion }
       coef:real;     { (i) the bandwidth expansion factor }
       length:integer      { (i) the length of lpc coefficient vectors }
   );
   procedure vq(
       Xq:pareal;      { (o) the quantized vector }
       index:pinteger;     { (o) the quantization index }
       CB:pareal;{ (i) the vector quantization codebook }
       X:pareal;       { (i) the vector to quantize }
       n_cb:integer;       { (i) the number of vectors in the codebook }
       dim:integer         { (i) the dimension of all vectors }
   );
   procedure SplitVQ(
       qX:PAreal;      { (o) the quantized vector }
       index:paInteger;     { (o) a vector of indexes for all vector
                              codebooks in the split }
       X:pareal;       { (i) the vector to quantize }
       CB:pareal;{ (i) the quantizer codebook }
       nsplit:integer;     { the number of vector splits }
       dim:PAInteger; { the dimension of X and qX }
       cbsize:PAInteger { the number of vectors in the codebook }
   );
   procedure sort_sq(
       xq:preal;      { (o) the quantized value }
       index:pinteger;     { (o) the quantization index }
       x:real;    { (i) the value to quantize }
       cb:pareal;{ (i) the quantization codebook }
       cb_size:integer      { (i) the size of the quantization codebook }
   );
   function LSF_check(    { (o) 1 for stable lsf vectors and 0 for
                              nonstable ones }
       lsf:pareal;     { (i) a table of lsf vectors }
       dim:integer;    { (i) the dimension of each lsf vector }
       NoAn:integer    { (i) the number of lsf vectors in the
                              table }
   ):integer;
implementation
   procedure autocorr(
       r:pareal;       { (o) autocorrelation vector }
       x:pareal; { (i) data vector }
       N:integer;          { (i) length of data vector }
       order:integer       { largest lag for calculated
                          autocorrelations }
   );
   var
       lag, nn:integer;
       sum:real;
   begin

       for lag := 0 to order do
       begin
           sum := 0;
           for nn := 0 to  N - lag-1 do
           begin
               sum :=sum+ x[nn] * x[nn+lag];
           end;
           r[lag] := sum;
       end;
   end;

   {----------------------------------------------------------------*
    *  window multiplication
    *---------------------------------------------------------------}

   procedure window(
       z:pareal;       { (o) the windowed data }
       x:pareal; { (i) the original data vector }
       y:pareal; { (i) the window }
       N:integer           { (i) length of all vectors }
   );
   var
   	i:integer;
   begin
       for i := 0 to N-1 do
       begin
           z[i] := x[i] * y[i];
       end;
   end;

   {----------------------------------------------------------------*
    *  levinson-durbin solution for lpc coefficients
    *---------------------------------------------------------------}

   procedure levdurb(
       a:pareal;       { (o) lpc coefficient vector starting
                              with 1.0 }
       k:pareal;       { (o) reflection coefficients }
       r:pareal;       { (i) autocorrelation vector }
       order:integer       { (i) order of lpc filter }
   );
   var
       sum, alpha:real;
       m, m_h, i:integer;
   begin

       a[0] := 1.0;

       if (r[0] < EPS) then 
       begin { if r[0] <:= 0, set LPC coeff. to zero }
           for i := 0 to  order-1 do
           begin
               k[i] := 0;
               a[i+1] := 0;
           end;
       end
       else 
       begin
           a[1] := -r[1]/r[0];
           k[0] := -r[1]/r[0];
           alpha := r[0] + r[1] * k[0];
           for m := 1 to order-1 do
           begin
               sum := r[m + 1];
               for i := 0 to m-1 do
               begin
                   sum :=sum + a[i+1] * r[m - i];
               end;
               k[m] := -sum / alpha;
               alpha :=alpha + k[m] * sum;
               m_h := (m + 1) shr 1;
               for i := 0 to m_h-1 do
               begin
                   sum := a[i+1] + k[m] * a[m - i];
                   a[m - i] :=a[m - i] + k[m] * a[i+1];
                   a[i+1] := sum;
               end;
               a[m+1] := k[m];
           end;
       end;
   end;

   {----------------------------------------------------------------*
    *  interpolation between vectors
    *---------------------------------------------------------------}

   procedure interpolate(
       nout:pareal;      { (o) the interpolated vector }
       in1:pareal;     { (i) the first vector for the
                              interpolation }
       in2:pareal;     { (i) the second vector for the
                              interpolation }
       coef:real;      { (i) interpolation weights }
       length:integer      { (i) length of all vectors }
   );
   var
       i:integer;
       invcoef:real;
   begin
       invcoef := 1.0 - coef;
       for i := 0 to length-1 do
       begin
           nout[i] := coef * in1[i] + invcoef * in2[i];
       end;
   end;

   {----------------------------------------------------------------*
    *  lpc bandwidth expansion
    *---------------------------------------------------------------}

   procedure bwexpand(
       nout:pareal;      { (o) the bandwidth expanded lpc
                              coefficients }
       nin:pareal;      { (i) the lpc coefficients before bandwidth
                              expansion }
       coef:real;     { (i) the bandwidth expansion factor }
       length:integer      { (i) the length of lpc coefficient vectors }
   );
   var
       i:integer;
       chirp:real;
   begin

       chirp := coef;

       nout[0] := nin[0];
       for i := 1 to length-1 do
       begin
           nout[i] := chirp * nin[i];
           chirp :=chirp * coef;
       end;
   end;

   {----------------------------------------------------------------*
    *  vector quantization
    *---------------------------------------------------------------}

   procedure vq(
       Xq:pareal;      { (o) the quantized vector }
       index:pinteger;     { (o) the quantization index }
       CB:pareal;{ (i) the vector quantization codebook }
       X:pareal;       { (i) the vector to quantize }
       n_cb:integer;       { (i) the number of vectors in the codebook }
       dim:integer         { (i) the dimension of all vectors }
   );
   var
       i, j:integer;
       pos, minindex:integer;
       dist, tmp, mindist:real;
   begin
       pos := 0;
       mindist := FLOAT_MAX;
       minindex := 0;
       for j := 0 to n_cb-1 do
       begin
           dist := X[0] - CB[pos];
           dist :=dist * dist;
           for i := 1 to dim-1 do
           begin
               tmp := X[i] - CB[pos + i];
               dist:=dist + tmp*tmp;
           end;

           if (dist < mindist) then
           begin
               mindist := dist;
               minindex := j;
           end;
           pos :=pos + dim;
       end;
       for i := 0 to dim-1 do
       begin
           Xq[i] := CB[minindex*dim + i];
       end;
       index^ := minindex;
   end;

   {----------------------------------------------------------------*
    *  split vector quantization
    *---------------------------------------------------------------}

   procedure SplitVQ(
       qX:PAreal;      { (o) the quantized vector }
       index:paInteger;     { (o) a vector of indexes for all vector
                              codebooks in the split }
       X:pareal;       { (i) the vector to quantize }
       CB:pareal;{ (i) the quantizer codebook }
       nsplit:integer;     { the number of vector splits }
       dim:PAInteger; { the dimension of X and qX }
       cbsize:PAInteger { the number of vectors in the codebook }
   );
   var
   	cb_pos, X_pos, i:integer;
   begin
       cb_pos := 0;
       X_pos:= 0;
       for i := 0 to nsplit-1 do
       begin
           vq(@qX[X_pos], @index[ i], @CB[ cb_pos], @X[X_pos],
               cbsize[i], dim[i]);
           X_pos :=X_pos+ dim[i];
           cb_pos :=cb_pos + dim[i] * cbsize[i];
       end;
   end;

   {----------------------------------------------------------------*
    *  scalar quantization
    *---------------------------------------------------------------}

   procedure sort_sq(
       xq:preal;      { (o) the quantized value }
       index:pinteger;     { (o) the quantization index }
       x:real;    { (i) the value to quantize }
       cb:pareal;{ (i) the quantization codebook }
       cb_size:integer      { (i) the size of the quantization codebook }
   );
   var
   	i:integer;
   begin
       if (x <= cb[0]) then
       begin
           index^ := 0;
           xq^ := cb[0];
       end
       else 
       begin
           i := 0;
           while ((x > cb[i]) and (i < cb_size - 1)) do
           begin
               inc(i);
           end;

           if (x > ((cb[i] + cb[i - 1])/2)) then
           begin
               index^ := i;
               xq^ := cb[i];
           end
           else
           begin
               index^ := i - 1;
               xq^ := cb[i - 1];
           end;
       end;
   end;

   {----------------------------------------------------------------*
    *  check for stability of lsf coefficients
    *---------------------------------------------------------------}

   function LSF_check(    { (o) 1 for stable lsf vectors and 0 for
                              nonstable ones }
       lsf:pareal;     { (i) a table of lsf vectors }
       dim:integer;    { (i) the dimension of each lsf vector }
       NoAn:integer    { (i) the number of lsf vectors in the
                              table }
   ):integer;
   var
     k,n,m, Nit, change,pos:integer;
     //tmp:real;
   begin
     Nit:=2;
     change:=0;

       { LSF separation check}

       for n:=0 to Nit-1 do 
       begin { Run through a couple of times }
           for m:=0 to NoAn-1 do
           begin { Number of analyses per frame }
               for k:=0 to (dim-2) do
               begin
                   pos:=m*dim+k;
                   if ((lsf[pos+1]-lsf[pos])<eps) then
                   begin
                       if (lsf[pos+1]<lsf[pos]) then
                       begin
                           //tmp:=lsf[pos+1];
                           lsf[pos+1]:= lsf[pos]+eps2;
                           lsf[pos]:= lsf[pos+1]-eps2;
                       end
                       else 
                       begin
                           lsf[pos]:=lsf[pos]-eps2;
                           lsf[pos+1]:=lsf[pos+1]-eps2;
                       end;
                       change:=1;
                   end;

                   if (lsf[pos]<minlsf) then
                   begin
                       lsf[pos]:=minlsf;
                       change:=1;
                   end;

                   if (lsf[pos]>maxlsf) then
                   begin
                       lsf[pos]:=maxlsf;
                       change:=1;
                   end;
               end;
           end;
       end;
       result:=change;
   end;
end.