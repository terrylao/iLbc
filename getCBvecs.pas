unit getCBvecs;

{$MODE Delphi}

interface
uses iLBC_define,constants,C2Delphi_header;

   {----------------------------------------------------------------*
    *  Construct codebook vector for given index.
    *---------------------------------------------------------------}
   procedure getCBvec(
       cbvec:pareal;   { (o) Constructed codebook vector }
       mem:pareal;     { (i) Codebook buffer }
       index:integer;      { (i) Codebook index }
       lMem:integer;       { (i) Length of codebook buffer }
       cbveclen:integer { (i) Codebook vector length }
   );
implementation
   procedure getCBvec(
       cbvec:pareal;   { (o) Constructed codebook vector }
       mem:pareal;     { (i) Codebook buffer }
       index:integer;      { (i) Codebook index }
       lMem:integer;       { (i) Length of codebook buffer }
       cbveclen:integer { (i) Codebook vector length }
   );
   var
       i,j, k, n, memInd, sFilt:integer;
       tmpbuf:array [0..CB_MEML-1] of real;
       base_size:integer;
       ilow, ihigh:integer;
       alfa, alfa1:real;
       tempbuff2:array [0..CB_MEML+CB_FILTERLEN] of real;
       pos:^real;
       pp, pp1:^real;
   begin

       { Determine size of codebook sections }

       base_size:=lMem-cbveclen+1;

       if (cbveclen=SUBL) then
       begin
           base_size:=base_size+cbveclen div 2;
       end;

       { No filter -> First codebook section }

       if (index<lMem-cbveclen+1) then
       begin
           { first non-interpolated vectors }
           k:=index+cbveclen;
           { get vector }
           move( mem[lMem-k],cbvec[0], cbveclen*sizeof(real));

       end
       else 
       if (index < base_size) then
       begin
           k:=2*(index-(lMem-cbveclen+1))+cbveclen;

           ihigh:=k div 2;
           ilow:=ihigh-5;

           { Copy first noninterpolated part }

           move( mem[lMem-k div 2],cbvec[0], ilow*sizeof(real));

           { interpolation }

           alfa1:=0.2;
           alfa:=0.0;
           for j:=ilow to ihigh-1 do
           begin
               cbvec[j]:=(1.0-alfa)*mem[lMem-k div 2+j]+alfa*mem[lMem-k+j];
               alfa:=alfa+alfa1;
           end;

           { Copy second noninterpolated part }

           move( mem[lMem-k+ihigh],cbvec[ihigh],(cbveclen-ihigh)*sizeof(real));

       end

       { Higher codebook section based on filtering }

       else 
       begin

           { first non-interpolated vectors }

           if (index-base_size<lMem-cbveclen+1) then
           begin

               fillchar(tempbuff2,CB_HALFFILTERLEN*sizeof(real), 0);
               move( mem[0],tempbuff2[CB_HALFFILTERLEN],lMem*sizeof(real));
               fillchar(tempbuff2[lMem+CB_HALFFILTERLEN],(CB_HALFFILTERLEN+1)*sizeof(real), 0);

               k:=index-base_size+cbveclen;
               sFilt:=lMem-k;
               memInd:=sFilt+1-CB_HALFFILTERLEN;

               { do filtering }
               pos:=@cbvec[0];
               fillchar(pos^, cbveclen*sizeof(real), 0);
               for n:=0 to cbveclen-1 do
               begin
                   pp:=@tempbuff2[memInd+n+CB_HALFFILTERLEN];
                   pp1:=@cbfiltersTbl[CB_FILTERLEN-1];
                   for j:=0 to CB_FILTERLEN-1 do
                   begin
                       pos^:=pos^+pp^*pp1^;
                       inc(pp);
                       dec(pp1);
                   end;
                   inc(pos);
               end;
           end

           { interpolated vectors }

           else
           begin
               fillchar(tempbuff2,CB_HALFFILTERLEN*sizeof(real), 0);
               move( mem[0],tempbuff2[CB_HALFFILTERLEN],lMem*sizeof(real));
               fillchar(tempbuff2[lMem+CB_HALFFILTERLEN],(CB_HALFFILTERLEN+1)*sizeof(real), 0);

               k:=2*(index-base_size-(lMem-cbveclen+1))+cbveclen;
               sFilt:=lMem-k;
               memInd:=sFilt+1-CB_HALFFILTERLEN;

               { do filtering }
               pos:=@tmpbuf[sFilt];
               fillchar(pos^, k*sizeof(real), 0);
               for i:=0 to k-1 do
               begin
                   pp:=@tempbuff2[memInd+i+CB_HALFFILTERLEN];
                   pp1:=@cbfiltersTbl[CB_FILTERLEN-1];
                   for j:=0 to CB_FILTERLEN-1 do 
                   begin
                       pos^:=pos^+pp^*pp1^;
                       inc(pp);
                       dec(pp1);
                   end;
                   inc(pos);
               end;

               ihigh:=k div 2;
               ilow:=ihigh-5;

               { Copy first noninterpolated part }

               move( tmpbuf[lMem-k div 2],cbvec[0],ilow*sizeof(real));

               { interpolation }

               alfa1:=0.2;
               alfa:=0.0;
               for j:=ilow to ihigh-1 do
               begin
                   cbvec[j]:=(1.0-alfa)*tmpbuf[lMem-k div 2+j]+alfa*tmpbuf[lMem-k+j];
                   alfa:=alfa+alfa1;
               end;

               { Copy second noninterpolated part }

               move( tmpbuf[lMem-k+ihigh],cbvec[ihigh],(cbveclen-ihigh)*sizeof(real));
           end;
       end;
   end;
end.