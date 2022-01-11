#pragma once
#include <cstdint>
#include <cstdio>
#include <exception>
/********************************************************************************
 * SP_dec: state classe for HT SigProp decoding
 *******************************************************************************/
class SP_dec {
 private:
  const uint32_t Lref;
  uint8_t bits;
  uint8_t tmp;
  uint8_t last;
  uint32_t pos;
  const uint8_t *Dref;

 public:
  SP_dec(const uint8_t *HT_magref_segment, uint32_t magref_length)
      : Lref(magref_length),
        bits(0),
        tmp(0),
        last(0),
        pos(0),
        Dref((Lref == 0) ? nullptr : HT_magref_segment) {}
  uint8_t importSigPropBit() {
    uint8_t val;
    if (bits == 0) {
      bits = (last == 0xFF) ? 7 : 8;
      if (pos < Lref) {
        tmp = *(Dref + pos);
        pos++;
        if ((tmp & (1 << bits)) != 0) {
          printf("ERROR: importSigPropBit error\n");
          throw std::exception();
        }
      } else {
        tmp = 0;
      }
      last = tmp;
    }
    val = tmp & 1;
    tmp >>= 1;
    bits--;
    return val;
  }
};

/********************************************************************************
 * MR_dec: state classe for HT MagRef decoding
 *******************************************************************************/
class MR_dec {
 private:
  const uint32_t Lref;
  uint8_t bits;
  uint8_t last;
  uint8_t tmp;
  int32_t pos;
  const uint8_t *Dref;

 public:
  MR_dec(const uint8_t *HT_magref_segment, uint32_t magref_length)
      : Lref(magref_length),
        bits(0),
        last(0xFF),
        tmp(0),
        pos((Lref == 0) ? -1 : magref_length - 1),
        Dref((Lref == 0) ? nullptr : HT_magref_segment) {}
  uint8_t importMagRefBit() {
    uint8_t val;
    if (bits == 0) {
      if (pos >= 0) {
        tmp = *(Dref + pos);
        pos--;
      } else {
        tmp = 0;
      }
      bits = 8;
      if (last > 0x8F && (tmp & 0x7F) == 0x7F) {
        bits = 7;
      }
      last = tmp;
    }
    val = tmp & 1;
    tmp >>= 1;
    bits--;
    return val;
  }
};