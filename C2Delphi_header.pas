unit C2Delphi_header;

{$MODE Delphi}

interface
const
	MaxArraySize=32767	;
type
	CharArrayType = array[0..MaxArraySize] of char;
	ByteArrayType = array[0..MaxArraySize] of byte;
	IntegerArrayType = array[0..MaxArraySize] of Integer;
	realArrayType = array[0..MaxArraySize] of real;
	CardinalArrayType = array[0..MaxArraySize] of Cardinal;
	shortArrayType = array[0..MaxArraySize] of Smallint;
	PAchar = ^CharArrayType;
	PAByte = ^ByteArrayType;
	PAInteger = ^IntegerArrayType;
	PAreal = ^realArrayType;
	PACardinal = ^CardinalArrayType;
  preal = ^real;
  PAshort = ^shortArrayType;
implementation
end.