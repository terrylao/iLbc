{
  publish with BSD Licence.
	Copyright (c) Terry Lao
}
unit lsf;

{$MODE Delphi}

interface
uses iLBC_define,C2Delphi_header;

   {----------------------------------------------------------------*
    *  conversion from lpc coefficients to lsf coefficients
    *---------------------------------------------------------------}

   procedure a2lsf(
       freq:pareal;{ (o) lsf coefficients }
       a:pareal    { (i) lpc coefficients }
   );

   {----------------------------------------------------------------*
    *  conversion from lsf coefficients to lpc coefficients
    *---------------------------------------------------------------}

   procedure lsf2a(
       a_coef:pareal;  { (o) lpc coefficients }
       freq:pareal     { (i) lsf coefficients }
   );


implementation
   procedure a2lsf(
       freq:pareal;{ (o) lsf coefficients }
       a:pareal    { (i) lpc coefficients }
   );
   var
       steps:array [0..LSF_NUMBER_OF_STEPS-1] of real;
       step:real;
       step_idx:integer;
       lsp_index:integer;
       p:array [0..LPC_HALFORDER-1] of real;
       q:array [0..LPC_HALFORDER-1] of real;
       p_pre:array [0..LPC_HALFORDER-1] of real;
       q_pre:array [0..LPC_HALFORDER-1] of real;
       old_p, old_q:real;
       old:^real;
       pq_coef:pareal;
       omega, old_omega:real;
       i:integer;
       hlp, hlp1, hlp2, hlp3, hlp4, hlp5:real;
   begin
			 steps[0]:=0.00635;
			 steps[1]:=0.003175;
			 steps[2]:=0.0015875;
			 steps[3]:=0.00079375;
       for i:=0 to LPC_HALFORDER-1 do
       begin
           p[i] := -1.0 * (a[i + 1] + a[LPC_FILTERORDER - i]);
           q[i] := a[LPC_FILTERORDER - i] - a[i + 1];
       end;

       p_pre[0] := -1.0 - p[0];
       p_pre[1] := - p_pre[0] - p[1];
       p_pre[2] := - p_pre[1] - p[2];
       p_pre[3] := - p_pre[2] - p[3];
       p_pre[4] := - p_pre[3] - p[4];
       p_pre[4] := p_pre[4] / 2;

       q_pre[0] := 1.0 - q[0];
       q_pre[1] := q_pre[0] - q[1];
       q_pre[2] := q_pre[1] - q[2];
       q_pre[3] := q_pre[2] - q[3];
       q_pre[4] := q_pre[3] - q[4];
       q_pre[4] := q_pre[4] / 2;

       omega := 0.0;
       old_omega := 0.0;

       old_p := FLOAT_MAX;
       old_q := FLOAT_MAX;

       { Here we loop through lsp_index to find all the
          LPC_FILTERORDER roots for omega. }

       for lsp_index := 0 to LPC_FILTERORDER-1 do
       begin
           { Depending on lsp_index being even or odd, we
           alternatively solve the roots for the two LSP equations. }

           if ((lsp_index and $1) = 0) then
           begin
               pq_coef := @p_pre;
               old := @old_p;
           end 
           else 
           begin
               pq_coef := @q_pre;
               old := @old_q;
           end;

           { Start with low resolution grid }
           step_idx := 0;
           step := steps[step_idx];
           while ( step_idx < LSF_NUMBER_OF_STEPS) do
           begin
               {  cos(10piw) + pq(0)cos(8piw) + pq(1)cos(6piw) +
               pq(2)cos(4piw) + pq(3)cod(2piw) + pq(4) }

               hlp := cos(omega * TWO_PI);
               hlp1 := 2.0 * hlp + pq_coef[0];
               hlp2 := 2.0 * hlp * hlp1 - 1.0 + pq_coef[1];
               hlp3 := 2.0 * hlp * hlp2 - hlp1 + pq_coef[2];
               hlp4 := 2.0 * hlp * hlp3 - hlp2 + pq_coef[3];
               hlp5 := hlp * hlp4 - hlp3 + pq_coef[4];

               if (((hlp5 * (old^)) <= 0.0) or (omega >= 0.5)) then
               begin
                   if (step_idx = (LSF_NUMBER_OF_STEPS - 1)) then
                   begin
                       if (abs(hlp5) >= abs(old^)) then
                       begin
                           freq[lsp_index] := omega - step;
                       end
                       else
                       begin
                           freq[lsp_index] := omega;
                       end;

                       if ((old^) >= 0.0) then
                       begin
                           old^ := -1.0 * FLOAT_MAX;
                       end
                       else
                       begin
                           old^ := FLOAT_MAX;
                       end;

                       omega := old_omega;
//                       step_idx := 0;

                       step_idx := LSF_NUMBER_OF_STEPS;
                   end
                   else
                   begin
                       if (step_idx = 0) then
                       begin
                           old_omega := omega;
                       end;

                       inc(step_idx);
                       omega :=omega - steps[step_idx];

                       { Go back one grid step }

                       step := steps[step_idx];
                   end;
               end
               else
               begin

               { increment omega until they are of different sign,
               and we know there is at least one root between omega
               and old_omega }
                   old^ := hlp5;
                   omega :=omega + step;
               end;
           end;
       end;

       for i := 0 to LPC_FILTERORDER-1 do
       begin
           freq[i] := freq[i] * TWO_PI;
       end;
   end;
   {----------------------------------------------------------------*
    *  conversion from lsf coefficients to lpc coefficients
    *---------------------------------------------------------------}

   procedure lsf2a(
       a_coef:pareal;  { (o) lpc coefficients }
       freq:pareal     { (i) lsf coefficients }
   );
   var
       i, j:integer;
       hlp:real;
       p:array [0..LPC_HALFORDER-1]  of real;
       q:array [0..LPC_HALFORDER] of real;
       a:array [0..LPC_HALFORDER] of real;
       a1:array [0..LPC_HALFORDER-1] of real;
       a2:array [0..LPC_HALFORDER-1] of real;
       b:array [0..LPC_HALFORDER] of real;
       b1:array [0..LPC_HALFORDER-1] of real;
       b2:array [0..LPC_HALFORDER-1] of real;
   begin
       for i:=0 to LPC_FILTERORDER-1 do
       begin
           freq[i] := freq[i] * PI2;
       end;

       { Check input for ill-conditioned cases.  This part is not
       found in the TIA standard.  It involves the following 2 IF
       blocks.  If "freq" is judged ill-conditioned, then we first
       modify freq[0] and freq[LPC_HALFORDER-1] (normally
       LPC_HALFORDER := 10 for LPC applications), then we adjust
       the other "freq" values slightly }


       if ((freq[0] <= 0.0) or (freq[LPC_FILTERORDER - 1] >= 0.5)) then
       begin
           if (freq[0] <= 0.0) then
           begin
               freq[0] := 0.022;
           end;

           if (freq[LPC_FILTERORDER - 1] >= 0.5) then
           begin
               freq[LPC_FILTERORDER - 1] := 0.499;
           end;

           hlp := (freq[LPC_FILTERORDER - 1] - freq[0]) / (LPC_FILTERORDER - 1);

           for i:=1 to LPC_FILTERORDER-1 do
           begin
               freq[i] := freq[i - 1] + hlp;
           end;
       end;

       fillchar(a1, LPC_HALFORDER*sizeof(real), 0);
       fillchar(a2, LPC_HALFORDER*sizeof(real), 0);
       fillchar(b1, LPC_HALFORDER*sizeof(real), 0);
       fillchar(b2, LPC_HALFORDER*sizeof(real), 0);
       fillchar(a , (LPC_HALFORDER+1)*sizeof(real), 0);
       fillchar(b , (LPC_HALFORDER+1)*sizeof(real), 0);

       { p[i] and q[i] compute cos(2*pi*omega_begin2jend;) and
       cos(2*pi*omega_begin2j-1end; in eqs. 4.2.2.2-1 and 4.2.2.2-2.
       Note that for this code p[i] specifies the coefficients
       used in .Q_A(z) while q[i] specifies the coefficients used
       in .P_A(z) }

       for i:=0 to LPC_HALFORDER-1 do
       begin
           p[i] := cos(TWO_PI * freq[2 * i]);
           q[i] := cos(TWO_PI * freq[2 * i + 1]);
       end;

       a[0] := 0.25;
       b[0] := 0.25;

       for i:= 0 to LPC_HALFORDER-1 do
       begin
           a[i + 1] := a[i] - 2 * p[i] * a1[i] + a2[i];
           b[i + 1] := b[i] - 2 * q[i] * b1[i] + b2[i];
           a2[i] := a1[i];
           a1[i] := a[i];
           b2[i] := b1[i];
           b1[i] := b[i];
       end;

       for j:=0 to LPC_FILTERORDER-1 do
       begin
           if (j = 0) then
           begin
               a[0] := 0.25;
               b[0] := -0.25;
           end
           else
           begin
               a[0] :=0.0;
               b[0] := 0.0;
           end;

           for i:=0 to LPC_HALFORDER-1 do
           begin
               a[i + 1] := a[i] - 2 * p[i] * a1[i] + a2[i];
               b[i + 1] := b[i] - 2 * q[i] * b1[i] + b2[i];
               a2[i] := a1[i];
               a1[i] := a[i];
               b2[i] := b1[i];
               b1[i] := b[i];
           end;
           a_coef[j + 1] := 2 * (a[LPC_HALFORDER] + b[LPC_HALFORDER]);
       end;
       a_coef[0] := 1.0;
   end;
end.
