program iLbc;
{$APPTYPE CONSOLE}
{$MODE Delphi}

uses
  SysUtils, Interfaces,
  iLBC_define in 'iLBC_define.pas',
  constants in 'constants.pas',
  C2Delphi_header in 'C2Delphi_header.pas',
  syntFilters in 'syntFilters.pas',
  lsf in 'lsf.pas',
  anaFilters in 'anaFilters.pas',
  createCB in 'createCB.pas',
  doCPLC in 'doCPLC.pas',
  filter in 'filter.pas',
  FrameClassifys in 'FrameClassifys.pas',
  getCBvecs in 'getCBvecs.pas',
  helpfun in 'helpfun.pas',
  hpInputs in 'hpInputs.pas',
  hpOutputs in 'hpOutputs.pas',
  gainquants in 'gainquants.pas',
  enhancers in 'enhancers.pas',
  StateSearchWs in 'StateSearchWs.pas',
  StateConstructWs in 'StateConstructWs.pas',
  iCBConstructs in 'iCBConstructs.pas',
  iCBSearchs in 'iCBSearchs.pas',
  LPCdecode in 'LPCdecode.pas',
  LPCencodes in 'LPCencodes.pas',
  packing in 'packing.pas',
  iLBC_decodes in 'iLBC_decodes.pas',
  iLBC_encodes in 'iLBC_encodes.pas';

{.$R *.res}

Const
   ILBCNOOFWORDS_MAX = (NO_OF_BYTES_30MS div 2);


   {---------------------------------------------------------------*
    *  Main program to test iLBC encoding and decoding
    *
    *  Usage:
    *    exefile_name.exe <infile> <bytefile> <outfile> <channel>
    *
    *    <infile>   : Input file, speech for encoder (16-bit pcm file)
    *    <bytefile> : Bit stream output from the encoder
    *    <outfile>  : Output file, decoded speech (16-bit pcm file)
    *    <channel>  : Bit error file, optional (16-bit)
    *                     1 - Packet received correctly
    *                     0 - Packet Lost
    *
    *--------------------------------------------------------------}
var
       starttime:real;
       runtime:real;
       outtime:real;

       ifileid,efileid,ofileid, cfileid,iBytesRead:integer;
       data:array [0..BLOCKL_MAX-1] of Smallint;
       encoded_data:array [0..ILBCNOOFWORDS_MAX-1] of Smallint;
       decoded_data:array [0..BLOCKL_MAX-1] of Smallint;
       len:integer;
       pli, mode:Smallint;
       blockcount:integer;
       packetlosscount:integer;
       { Create structs }
       Enc_Inst:iLBC_Enc_Inst_t;
       Dec_Inst:iLBC_Dec_Inst_t;
   {----------------------------------------------------------------*
    *  Encoder interface function
    *---------------------------------------------------------------}

   function encode(   { (o) Number of bytes encoded }
       iLBCenc_inst:piLBC_Enc_Inst_t;
                                   { (i/o) Encoder instance }
       encoded_data:PAshort;    { (o) The encoded bytes }
       data:PAshort                 { (i) The signal block to encode}
   ):Smallint;
   var
       block:array [0..BLOCKL_MAX-1] of real;
       k:integer;
   begin

       { convert signal to float }

       for k:=0 to iLBCenc_inst^.blockl-1 do
           block[k] := data[k];

       { do the actual encoding }

       iLBC_encode(@encoded_data[0], @block, iLBCenc_inst);

       result:= Smallint(iLBCenc_inst^.no_of_bytes);
   end;

   {----------------------------------------------------------------*
    *  Decoder interface function
    *---------------------------------------------------------------}

   function decode(       { (o) Number of decoded samples }
       iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) Decoder instance }
       decoded_data:PAshort;        { (o) Decoded signal block}
       encoded_data:PAshort;        { (i) Encoded bytes }
       mode:Smallint                       { (i) 0:=PL, 1:=Normal }
   ):Smallint;
   var
       k:integer;
       decblock:array [0..BLOCKL_MAX-1] of real;
       dtmp:real;
   begin

       { check if mode is valid }

       if (mode<0)  or  (mode>1) then
       begin
           writeln(ErrOutput,'ERROR - Wrong mode - 0, 1 allowed');
           result:=3;
           exit;
       end;

       { do actual decoding of block }

       iLBC_decode(@decblock, @encoded_data[0], iLBCdec_inst, mode);

       { convert to short }

       for k:=0 to iLBCdec_inst^.blockl-1 do
       begin
           dtmp:=decblock[k];

           if (dtmp<MIN_SAMPLE) then
               dtmp:=MIN_SAMPLE
           else
           if (dtmp>MAX_SAMPLE) then
               dtmp:=MAX_SAMPLE;
           decoded_data[k] := trunc(dtmp);
       end;

       result:= Smallint(iLBCdec_inst^.blockl);
   end;
   begin
   try
       { Runtime statistics }

       blockcount := 0;
       packetlosscount := 0;

       { get arguments and open files }

       if ((ParamCount<>4)  and  (ParamCount<>5)) then
       begin
           writeln(ErrOutput,inttostr(ParamCount));
           writeln(ErrOutput,
           '*-----------------------------------------------*');
           writeln(ErrOutput,Format('   %s <20,30> input encoded decoded (channel)',[ParamStr(0)]));
           writeln(ErrOutput,
           '   mode    : Frame size for the encoding/decoding');
           writeln(ErrOutput,
           '                 20 - 20 ms');
           writeln(ErrOutput,
           '                 30 - 30 ms');
           writeln(ErrOutput,
           '   input   : Speech for encoder (16-bit pcm file)');
           writeln(ErrOutput,
           '   encoded : Encoded bit stream');
           writeln(ErrOutput,
           '   decoded : Decoded speech (16-bit pcm file)');
           writeln(ErrOutput,
           '   channel : Packet loss pattern, optional (16-bit)');
           writeln(ErrOutput,
           '                  1 - Packet received correctly');
           writeln(ErrOutput,
           '                  0 - Packet Lost');
           writeln(ErrOutput,
           '*-----------------------------------------------*');
           exit;
       end;
       mode:=strtoint(ParamStr(1));
       if (mode <> 20)  and  (mode <> 30) then
       begin
           writeln(ErrOutput,Format('Wrong mode %s, must be 20, or 30',[ParamStr(1)]));
           exit;
       end;
       ifileid := FileOpen(ParamStr(2), fmOpenRead or fmShareDenyNone);

       if ( ifileid < 1) then
       begin
           writeln(ErrOutput,Format('Cannot open input file %s', [ParamStr(2)]));
           exit;
       end;
       if FileExists(ParamStr(3)) then
       begin
        SysUtils.DeleteFile(ParamStr(3));
       end;
       efileid := FileCreate(ParamStr(3));
       if ( efileid<1) then
       begin
           writeln(ErrOutput,Format('Cannot open encoded file %s',
               [ParamStr(3)]));
               exit;
       end;
       if FileExists(ParamStr(4)) then
       begin
        SysUtils.DeleteFile(ParamStr(4));
       end;
       ofileid := FileCreate(ParamStr(4));
       if (ofileid<1) then
       begin
           writeln(ErrOutput,Format('Cannot open decoded file %s',
               [ParamStr(4)]));
           exit;
       end;
       if (ParamCount=5) then
       begin
       		cfileid := FileOpen(ParamStr(5), fmOpenRead or fmShareDenyNone);
           if(cfileid<1) then
           begin
               writeln(ErrOutput,Format('Cannot open channel file %s',
                   [ParamStr(5)]));
               exit;
           end;
       end
       else
       begin
           cfileid:=0;
       end;

       { print info }

       writeln(ErrOutput,'');
       writeln(ErrOutput,
           '*---------------------------------------------------*');
       writeln(ErrOutput,
           '*                                                   *');
       writeln(ErrOutput,
           '*      iLBC test program                            *');
       writeln(ErrOutput,
           '*                                                   *');
       writeln(ErrOutput,
           '*                                                   *');
       writeln(ErrOutput,
           '*---------------------------------------------------*');
       writeln(ErrOutput,Format('Mode           : %2d ms', [mode]));
       writeln(ErrOutput,Format('Input file     : %s', [ParamStr(2)]));
       writeln(ErrOutput,Format('Encoded file   : %s', [ParamStr(3)]));
       writeln(ErrOutput,Format('Output file    : %s', [ParamStr(4)]));
       if (ParamCount =5) then
       begin
           writeln(ErrOutput,Format('Channel file   : %s', [ParamStr(5)]));
       end;
       writeln(ErrOutput,'');

       { Initialization }

       initEncode(@Enc_Inst, mode);
       initDecode(@Dec_Inst, mode, 1);

       { Runtime statistics }

       starttime:=GetTickCount();

       { loop over input blocks }
       iBytesRead := FileRead(ifileid, data[0], sizeof(Smallint)*Enc_Inst.blockl);
       while (iBytesRead=sizeof(Smallint)*Enc_Inst.blockl) do
       begin

           inc(blockcount);

           { encoding }

           write(ErrOutput,Format('--- Encoding block %d --- ',[blockcount]));
           len:=encode(@Enc_Inst, @encoded_data, @data);
           write(ErrOutput,#13);

           { write byte file }
           FileWrite(efileid, encoded_data[0], sizeof(char)*len);

           { get channel data if provided }
           if (ParamCount =6) then
           begin
               if fileread(cfileid,pli, sizeof(Smallint))<1 then
               begin
                   if ((pli<>0) and (pli<>1)) then
                   begin
                       writeln(ErrOutput,'Error in channel file');
                       exit;
                   end;
                   if (pli=0) then
                   begin
                       { Packet loss ^. remove info from frame }
                       fillchar(encoded_data,sizeof(Smallint)*ILBCNOOFWORDS_MAX, 0);
                       inc(packetlosscount);
                   end;
               end
               else
               begin
                   writeln(ErrOutput,'Error. Channel file too short');
                   exit;
               end;
           end
           else
           begin
               pli:=1;
           end;

           { decoding }

           write(ErrOutput,Format('--- Decoding block %d --- ',[blockcount]));

           len:=decode(@Dec_Inst, @decoded_data, @encoded_data, pli);
           write(ErrOutput,#13);

           { write output file }

           filewrite(ofileid,decoded_data,sizeof(Smallint)*len);
           iBytesRead := FileRead(ifileid, data[0], sizeof(Smallint)*Enc_Inst.blockl);
       end;

       { Runtime statistics }

       runtime := (GetTickCount()-starttime);
       outtime := (blockcount*mode/1000.0);
       writeln(ErrOutput,Format('Length of speech file: %.1f s', [outtime]));
       writeln(ErrOutput,Format('Packet loss          : %.1f%%',
           [100.0*packetlosscount/blockcount]));

       write('Time to run iLBC     :');
       writeln(ErrOutput,Format(' %.1f s (%.1f %% of realtime)', [runtime/1000,(100*runtime/outtime)/1000]));

       { close files }
 except
       FileClose(ifileid);  FileClose(efileid); FileClose(ofileid);
       if (ParamCount=6) then
       begin
           FileClose(cfileid);
       end;
 end;
end.
