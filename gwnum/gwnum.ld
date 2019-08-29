SECTIONS
{
  .data : { . = ALIGN(0x20); *(_GWDATA) }
}
