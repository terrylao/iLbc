unit iLbcCodec;
interface
uses
  C2Delphi_header,
  SysUtils,windows,
  anaFilters,
  constants,
  createCB,
  doCPLC,
  enhancers,
  filter,
  FrameClassifys,
  gainquants,
  getCBvecs,
  helpfun,
  hpInputs,
  hpOutputs,
  iCBConstructs,
  iCBSearchs,
  iLBC_decodes,
  iLBC_define,
  iLBC_encodes;

Const
  ILBCNOOFWORDS_MAX = (NO_OF_BYTES_30MS div 2);
	BITRATE_30MS	= NO_OF_BYTES_30MS*8*8000/BLOCKL_30MS; //raw bits per second = BytesPerFrame * BitsPerSample*SamplesPerSec/SamplesPerFrame
	BITRATE_20MS	= NO_OF_BYTES_20MS*8*8000/BLOCKL_20MS;
var
       data:array [0..BLOCKL_MAX-1] of Smallint;
       encoded_data:array [0..ILBCNOOFWORDS_MAX-1] of Smallint;
       decoded_data:array [0..BLOCKL_MAX-1] of Smallint;
       pli, mode:Smallint;
       { Create structs }
       Enc_Inst:iLBC_Enc_Inst_t;
       Dec_Inst:iLBC_Dec_Inst_t;
function encode(   { (o) Number of bytes encoded }
   iLBCenc_inst:piLBC_Enc_Inst_t;
                               { (i/o) Encoder instance }
   encoded_data:PAshort;    { (o) The encoded bytes }
   encoded_DataByteSize:cardinal;
   data:PAshort;                 { (i) The signal block to encode}
   dataByteCount:cardinal
):cardinal;
function decode(       { (o) Number of decoded samples }
   iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) Decoder instance }
   decoded_data:PAshort;        { (o) Decoded signal block}
   decodedBufByteSize:cardinal;
   encoded_data:PAshort;        { (i) Encoded bytes }
   encodedByteCount:cardinal;
   mode:Smallint                       { (i) 0:=PL, 1:=Normal }
):cardinal;
   function encode_single(   { (o) Number of bytes encoded }
       iLBCenc_inst:piLBC_Enc_Inst_t;
                                   { (i/o) Encoder instance }
       encoded_data:PAshort;    { (o) The encoded bytes }
       data:PAshort                 { (i) The signal block to encode}
   ):Smallint;
   function decode_single(       { (o) Number of decoded samples }
       iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) Decoder instance }
       decoded_data:PAshort;        { (o) Decoded signal block}
       encoded_data:PAshort;        { (i) Encoded bytes }
       mode:Smallint                       { (i) 0:=PL, 1:=Normal }
   ):Smallint;
   
implementation
   function encode_single(   { (o) Number of bytes encoded }
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

   function decode_single(       { (o) Number of decoded samples }
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
   {----------------------------------------------------------------*
    *  Encoder interface function
    *---------------------------------------------------------------}

function encode(   { (o) Number of bytes encoded }
   iLBCenc_inst:piLBC_Enc_Inst_t;
                               { (i/o) Encoder instance }
   encoded_data:PAshort;    { (o) The encoded bytes }
   encoded_DataByteSize:cardinal;
   data:PAshort;                 { (i) The signal block to encode}
   dataByteCount:cardinal
):cardinal;{肚琌byte 计, 传Θshort int 璶埃 2}
var
   block:array [0..BLOCKL_MAX-1] of real;
   k,i,j:integer;
begin//240 short int(480 bytes) 溃Θ25 short int(50 byte)
	j:=dataByteCount div sizeof(short);
	j:=j - (j mod iLBCenc_inst^.blockl);
	i:=0;
	result:=0;
	while i<j do
	begin
		{ convert signal to float }
		
		for k:=0 to iLBCenc_inst^.blockl-1 do
		   block[k] := data[k+i];
		i:=i+iLBCenc_inst^.blockl;
		{ do the actual encoding }
		
		iLBC_encode(@encoded_data[result div 2], @block, iLBCenc_inst);
		
		result:=result + Smallint(iLBCenc_inst^.no_of_bytes);
	end;
	j:=(dataByteCount div sizeof(short)) mod iLBCenc_inst^.blockl;
	if j>0 then
	begin
		fillchar(block[0],BLOCKL_MAX*sizeof(real),0);
		for k:=0 to j-1 do
		   block[k] := data[k+i];
		
		iLBC_encode(@encoded_data[result div 2], @block, iLBCenc_inst);
		
		result:=result + Smallint(iLBCenc_inst^.no_of_bytes);
	end;
end;

   {----------------------------------------------------------------*
    *  Decoder interface function
    *---------------------------------------------------------------}

function decode(       { (o) Number of decoded samples }
   iLBCdec_inst:piLBC_Dec_Inst_t;  { (i/o) Decoder instance }
   decoded_data:PAshort;        { (o) Decoded signal block}
   decodedBufByteSize:cardinal;
   encoded_data:PAshort;        { (i) Encoded bytes }
   encodedByteCount:cardinal;
   mode:Smallint                       { (i) 0:=PL, 1:=Normal }
):cardinal;{return 琌short int 计ヘ, 传Θbyte 璶 2}
var
   k,i,j:integer;
   decblock:array [0..BLOCKL_MAX-1] of real;
   dtmp:real;
begin

	{ check if mode is valid }
	
	if (mode<0)  or  (mode>1) then
	begin
		result:=0;
		exit;
	end;
	//50 byte 25 bytes --> 240 short int480 bytes
	{ do actual decoding of block }
	j:=encodedByteCount div 2;
	i:=0;
  result:=0;
	while i<j do
	begin
		iLBC_decode(@decblock, @encoded_data[i], iLBCdec_inst, mode);

		{ convert to short }

		for k:=0 to iLBCdec_inst^.blockl-1 do
		begin
			dtmp:=decblock[k];

			if (dtmp<MIN_SAMPLE) then
			   dtmp:=MIN_SAMPLE
			else
			if (dtmp>MAX_SAMPLE) then
			   dtmp:=MAX_SAMPLE;
			decoded_data[k+result] := trunc(dtmp);
		end;

		inc(i,25);
		result:=result + cardinal(iLBCdec_inst^.blockl);
	end;
	j:=encodedByteCount mod 50;
	if j>0 then
	begin
		move(encoded_data[i],encoded_data[i-((50-j) div 2)],j);
		fillchar(encoded_data[i+(j div 2)],(50-j) div 2,0);
		iLBC_decode(@decblock, @encoded_data[i-((50-j) div 2)], iLBCdec_inst, mode);
		
		{ convert to short }
		
		for k:=0 to iLBCdec_inst^.blockl-1 do
		begin
			dtmp:=decblock[k];
			
			if (dtmp<MIN_SAMPLE) then
			   dtmp:=MIN_SAMPLE
			else
			if (dtmp>MAX_SAMPLE) then
			   dtmp:=MAX_SAMPLE;
			decoded_data[k+result] := trunc(dtmp);
		end;

		result:=result + cardinal(iLBCdec_inst^.blockl);	
	end;
end;
begin
  initEncode(@Enc_Inst, 30);
  initDecode(@Dec_Inst, 30, 1);
end.

 
