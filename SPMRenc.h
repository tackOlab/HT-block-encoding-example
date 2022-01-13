#pragma once
#include <cstdint>
#include <cstdio>
#include <exception>
class MR_enc;  // forward declaration for friend function termSPandMR()
/********************************************************************************
 * SP_enc: state classe for HT SigProp encoding
 *******************************************************************************/
class SP_enc {
 private:
  uint32_t pos;
  uint8_t bits;
  uint8_t max;
  uint8_t tmp;
  uint8_t *const buf;
  friend int32_t termSPandMR(SP_enc &, MR_enc &);

 public:
  SP_enc(uint8_t *Dref) : pos(0), bits(0), max(8), tmp(0), buf(Dref) {}
  void emitSPBit(uint8_t bit) {
    tmp |= (bit << bits);
    bits++;
    if (bits == max) {
      buf[pos] = tmp;
      pos++;
      max  = (tmp == 0xFF) ? 7 : 8;
      tmp  = 0;
      bits = 0;
    }
  }
  void termSP() {
    if (tmp != 0) {
      buf[pos] = tmp;
      pos++;
      max = (tmp == 0xFF) ? 7 : 8;
    }
    if (max == 7) {
      buf[pos] = 0x00;
      pos++;  // this prevents the appearance of a terminal 0xFF
    }
  }
  uint32_t getPos() { return pos; }
  void show() {
    for (int i = 0; i < pos; ++i) {
      printf("%02x ", buf[i]);
    }
    printf("\n");
  }
};
/********************************************************************************
 * MR_enc: state classe for HT MagRef encoding
 *******************************************************************************/
class MR_enc {
 private:
  uint32_t pos;
  uint8_t bits;
  uint8_t tmp;
  uint8_t last;
  uint8_t *const buf;
  friend int32_t termSPandMR(SP_enc &, MR_enc &);

 public:
  MR_enc(uint8_t *Dref) : pos(0), bits(0), tmp(0), last(255), buf(Dref) {}
  void emitMRBit(uint8_t bit) {
    tmp |= (bit << bits);
    bits++;
    if ((last > 0x8F) && (tmp == 0x7F)) {
      bits++;  // this must leave MR_bits equal to 8
    }
    if (bits == 8) {
      buf[pos] = tmp;
      pos++;
      last = tmp;
      tmp  = 0;
      bits = 0;
    }
  }
  uint32_t getPos() { return pos; }
  void show() {
    for (int i = 0; i < pos; ++i) {
      printf("%02x ", buf[i]);
    }
    printf("\n");
  }
};

int32_t termSPandMR(SP_enc &SP, MR_enc &MR) {
  uint8_t SP_mask = 0xFF >> (8 - SP.bits);  // if SP_bits is 0, SP_mask = 0
  SP_mask |= ((1 << SP.max) & 0x80);        // Auguments SP_mask to cover any stuff bit
  uint8_t MR_mask = 0xFF >> (8 - MR.bits);  // if MR_bits is 0, MR_mask = 0
  if ((SP_mask | MR_mask) == 0) {
    return 0;  // last SP byte cannot be 0xFF, since then SP_max would be 7
  }
  uint8_t fuse = SP.tmp | MR.tmp;
  if ((((fuse ^ SP.tmp) & SP_mask) | ((fuse ^ MR.tmp) & MR_mask)) == 0) {
    SP.buf[SP.pos] = fuse;  // fuse always < 0x80 here; no false marker risk
  } else {
    SP.buf[SP.pos] = SP.tmp;  // SP_tmp cannot be 0xFF
    MR.buf[MR.pos] = MR.tmp;
    MR.pos++;
  }
  SP.pos++;
  for (int i = MR.pos - 1, j = 0; i >= 0; --i, ++j) {
    SP.buf[SP.pos + j] = MR.buf[i];
  }
  return 1;
}
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