#include "SPMRenc.h"
#include "codingunits.h"

auto process_stripes_block = [](SP_enc &SigProp, j2k_codeblock *block, const uint16_t i_start,
                                const uint16_t j_start, const uint16_t width, const uint16_t height,
                                const uint16_t dum_stride) {
  uint8_t *sp;
  uint8_t causal_cond = 0;
  uint8_t bit;
  uint8_t mbr;
  uint32_t mbr_info;  // NOT USED

  for (int16_t j = j_start; j < j_start + width; j++) {
    mbr_info = 0;
    for (int16_t i = i_start; i < i_start + height; i++) {
      sp          = &block->block_states[(j + 1) + (i + 1) * (block->size.x + 2)];
      causal_cond = (((block->Cmodes & CAUSAL) == 0) || (i != i_start + height - 1));
      mbr         = 0;
      if (block->get_state(Sigma, i, j) == 0) {
        block->calc_mbr(mbr, i, j, mbr_info & 0x1EF, causal_cond);
      }
      mbr_info >>= 3;
      if (mbr != 0) {
        bit = (*sp >> 5) & 1;
        SigProp.emitSPBit(bit);
        block->modify_state(refinement_indicator, 1, i, j);
        block->modify_state(refinement_value, bit, i, j);
      }
      block->modify_state(scan, 1, i, j);
    }
  }
  for (uint16_t j = j_start; j < j_start + width; j++) {
    for (uint16_t i = i_start; i < i_start + height; i++) {
      sp = &block->block_states[(j + 1) + (i + 1) * (block->size.x + 2)];
      // encode sign
      if (block->get_state(Refinement_value, i, j)) {
        bit = (*sp >> 6) & 1;
        SigProp.emitSPBit(bit);
      }
    }
  }
};

void ht_sigprop_encode(j2k_codeblock *block, SP_enc &SigProp) {
  const uint16_t num_v_stripe = block->size.y / 4;
  const uint16_t num_h_stripe = block->size.x / 4;
  uint16_t i_start            = 0, j_start;
  uint16_t width              = 4;
  uint16_t width_last;
  uint16_t height           = 4;
  const uint16_t dum_stride = block->size.x + 2;

  // encode full-height (=4) stripes
  for (uint16_t n1 = 0; n1 < num_v_stripe; n1++) {
    j_start = 0;
    for (uint16_t n2 = 0; n2 < num_h_stripe; n2++) {
      process_stripes_block(SigProp, block, i_start, j_start, width, height, dum_stride);
      j_start += 4;
    }
    width_last = block->size.x % 4;
    if (width_last) {
      process_stripes_block(SigProp, block, i_start, j_start, width_last, height, dum_stride);
    }
    i_start += 4;
  }
  // encode remaining height stripes
  height  = block->size.y % 4;
  j_start = 0;
  for (uint16_t n2 = 0; n2 < num_h_stripe; n2++) {
    process_stripes_block(SigProp, block, i_start, j_start, width, height, dum_stride);
    j_start += 4;
  }
  width_last = block->size.x % 4;
  if (width_last) {
    process_stripes_block(SigProp, block, i_start, j_start, width_last, height, dum_stride);
  }
}

void ht_magref_encode(j2k_codeblock *block, MR_enc &MagRef) {
  const uint16_t blk_height   = block->size.y;
  const uint16_t blk_width    = block->size.x;
  const uint16_t num_v_stripe = block->size.y / 4;
  uint16_t i_start            = 0;
  uint16_t height             = 4;
  uint8_t *sp;
  uint8_t bit;

  for (uint16_t n1 = 0; n1 < num_v_stripe; n1++) {
    for (uint16_t j = 0; j < blk_width; j++) {
      for (uint16_t i = i_start; i < i_start + height; i++) {
        sp = &block->block_states[(j + 1) + (i + 1) * (block->size.x + 2)];
        if (block->get_state(Sigma, i, j) != 0) {
          bit = (sp[0] >> 5) & 1;
          MagRef.emitMRBit(bit);
          block->modify_state(refinement_indicator, 1, i, j);
        }
      }
    }
    i_start += 4;
  }
  height = blk_height % 4;
  for (uint16_t j = 0; j < blk_width; j++) {
    for (uint16_t i = i_start; i < i_start + height; i++) {
      sp = &block->block_states[(j + 1) + (i + 1) * (block->size.x + 2)];
      if (block->get_state(Sigma, i, j) != 0) {
        bit = (sp[0] >> 5) & 1;
        MagRef.emitMRBit(bit);
        block->modify_state(refinement_indicator, 1, i, j);
      }
    }
  }
}

int32_t termSPandMR(SP_enc &SP, MR_enc &MR) {
  uint8_t SP_mask = 0xFF >> (8 - SP.bits);  // if SP_bits is 0, SP_mask = 0
  SP_mask |= ((1 << SP.max) & 0x80);        // Auguments SP_mask to cover any stuff bit
  uint8_t MR_mask = 0xFF >> (8 - MR.bits);  // if MR_bits is 0, MR_mask = 0
  if ((SP_mask | MR_mask) == 0) {
    // last SP byte cannot be 0xFF, since then SP_max would be 7
    memmove(&SP.buf[SP.pos], &MR.buf[MR.pos + 1], MAX_Lref - MR.pos);
    return SP.pos + MAX_Lref - MR.pos;
  }
  uint8_t fuse = SP.tmp | MR.tmp;
  if ((((fuse ^ SP.tmp) & SP_mask) | ((fuse ^ MR.tmp) & MR_mask)) == 0) {
    SP.buf[SP.pos] = fuse;  // fuse always < 0x80 here; no false marker risk
  } else {
    SP.buf[SP.pos] = SP.tmp;  // SP_tmp cannot be 0xFF
    MR.buf[MR.pos] = MR.tmp;
    MR.pos--;  // MR buf gorws reverse order
  }
  SP.pos++;
  // MRのデータは左右反転してSPの後ろにくっつく
  memmove(&SP.buf[SP.pos], &MR.buf[MR.pos + 1], MAX_Lref - MR.pos);
  return SP.pos + MAX_Lref - MR.pos;
}