#pragma once
#include <cstdint>

#define MAX_Lcup 16834
#define MAX_Scup 4079
#define MAX_Lref 2046

/********************************************************************************
 * state_MS_enc: state class for MagSgn encoding
 *******************************************************************************/
class state_MS_enc {
 private:
  int32_t pos;         // current position in the buffer
  uint8_t *const buf;  // buffer for MagSgn
#ifdef MSNAIVE
  uint8_t bits;
  uint8_t max;
  uint8_t tmp;
#else
  uint64_t Creg;      // temporal buffer to store up to 4 codewords
  uint32_t ctreg;     // number of used bits in Creg
  uint8_t last;       // last byte in the buffer
  void emit_dword();  // internal function to emit 4 code words
#endif

 public:
  explicit state_MS_enc(uint8_t *p)
      : pos(0),
        buf(p),
#ifdef MSNAIVE
        bits(0),
        max(8),
        tmp(0)
#else
        Creg(0),
        ctreg(0),
        last(0)
#endif
  {
  }
#ifdef MSNAIVE
  void emitMagSgnBits(uint32_t cwd, uint8_t m_n);
#else
  void emitMagSgnBits(uint32_t cwd, uint8_t m_n, uint8_t emb_1);
#endif
  int32_t termMS();
};

class state_MEL_enc;  // forward declaration for friend function "termMELandVLC()"
/********************************************************************************
 * state_VLC_enc: state class for VLC encoding
 *******************************************************************************/
class state_VLC_enc {
 private:
  int32_t pos;
  uint8_t bits;
  uint8_t tmp;
  uint8_t last;
  uint8_t *const buf;

  friend int32_t termMELandVLC(state_VLC_enc &, state_MEL_enc &);

 public:
  explicit state_VLC_enc(uint8_t *p) : pos(MAX_Scup - 2), bits(4), tmp(0xF), last(0xFF), buf(p) {
    buf[pos + 1] = 0xFF;
  }
  void emitVLCBits(uint16_t cwd, uint8_t len);
};

/********************************************************************************
 * state_MEL_enc: state class for MEL encoding
 *******************************************************************************/
class state_MEL_enc {
 private:
  int8_t MEL_k;
  uint8_t MEL_run;
  const uint8_t MEL_E[13];
  uint8_t MEL_t;
  int32_t pos;
  uint8_t rem;
  uint8_t tmp;
  uint8_t *const buf;
  void emitMELbit(uint8_t bit);

  friend int32_t termMELandVLC(state_VLC_enc &, state_MEL_enc &);

 public:
  explicit state_MEL_enc(uint8_t *p)
      : MEL_k(0),
        MEL_run(0),
        MEL_E{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5},
        MEL_t(static_cast<uint8_t>(1 << MEL_E[MEL_k])),
        pos(0),
        rem(8),
        tmp(0),
        buf(p) {}
  void encodeMEL(uint8_t smel);
  void termMEL();
};