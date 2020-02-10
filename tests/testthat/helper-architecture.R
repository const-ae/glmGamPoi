
is_32_bit_architecture <- function(){
  .Machine$sizeof.pointer == 4
}
