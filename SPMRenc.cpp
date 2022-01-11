#include "SPMRenc.h"

#include "codingunits.h"

void j2k_codeblock::set_MagSgn_and_sigma(uint32_t &or_val) {
  const uint32_t height = this->size.y;
  const uint32_t width = this->size.x;
  const uint32_t stride = this->band_stride;

  for (uint16_t i = 0; i < height; ++i) {
    sprec_t *const sp = this->i_samples + i * stride;
    int32_t *const dp = this->sample_buf.get() + i * width;
    size_t block_index = (i + 1) * (size.x + 2) + 1;
    for (uint16_t j = 0; j < width; ++j) {
      int32_t temp = sp[j];
      uint32_t sign = static_cast<uint32_t>(temp) & 0x80000000;
      if (temp) {
        or_val |= 1;
        block_states[block_index] |= 1;
        // convert sample value to MagSgn
        temp = (temp < 0) ? -temp : temp;
        temp &= 0x7FFFFFFF;
        temp--;
        temp <<= 1;
        temp += sign >> 31;
        dp[j] = temp;
      }
      block_index++;
    }
  }
}

auto process_stripes_block =
    [](SP_dec &SigProp, j2k_codeblock *block, const uint16_t i_start,
       const uint16_t j_start, const uint16_t width, const uint16_t height,
       const uint16_t dum_stride, const uint8_t &pLSB) {
      int32_t *sp;
      uint8_t causal_cond = 0;
      uint8_t bit;
      uint8_t mbr;
      uint32_t mbr_info;  // NOT USED

      for (int16_t j = j_start; j < j_start + width; j++) {
        mbr_info = 0;
        for (int16_t i = i_start; i < i_start + height; i++) {
          sp = &block->sample_buf[j + i * block->size.x];
          causal_cond =
              (((block->Cmodes & CAUSAL) == 0) || (i != i_start + height - 1));
          mbr = 0;
          if (block->get_state(Sigma, i, j) == 0) {
            block->calc_mbr(mbr, i, j, mbr_info & 0x1EF, causal_cond);
          }
          mbr_info >>= 3;
          if (mbr != 0) {
            block->modify_state(refinement_indicator, 1, i, j);
            bit = SigProp.importSigPropBit();
            block->modify_state(refinement_value, bit, i, j);
            // block->set_refinement_value(bit, i, j);
            *sp |= bit << pLSB;
          }
          block->modify_state(scan, 1, i, j);
          // block->update_scan_state(1, i, j);
        }
      }
      for (uint16_t j = j_start; j < j_start + width; j++) {
        for (uint16_t i = i_start; i < i_start + height; i++) {
          sp = &block->sample_buf[j + i * block->size.x];
          // decode sign
          if ((*sp & (1 << pLSB)) != 0) {
            *sp = (*sp & 0x7FFFFFFF) | (SigProp.importSigPropBit() << 31);
          }
        }
      }
    };

void ht_sigprop_decode(j2k_codeblock *block, uint8_t *HT_magref_segment,
                       uint32_t magref_length, const uint8_t &pLSB) {
  SP_dec SigProp(HT_magref_segment, magref_length);
  const uint16_t num_v_stripe = block->size.y / 4;
  const uint16_t num_h_stripe = block->size.x / 4;
  uint16_t i_start = 0, j_start;
  uint16_t width = 4;
  uint16_t width_last;
  uint16_t height = 4;
  const uint16_t dum_stride = block->size.x + 2;

  // decode full-height (=4) stripes
  for (uint16_t n1 = 0; n1 < num_v_stripe; n1++) {
    j_start = 0;
    for (uint16_t n2 = 0; n2 < num_h_stripe; n2++) {
      process_stripes_block(SigProp, block, i_start, j_start, width, height,
                            dum_stride, pLSB);
      j_start += 4;
    }
    width_last = block->size.x % 4;
    if (width_last) {
      process_stripes_block(SigProp, block, i_start, j_start, width_last,
                            height, dum_stride, pLSB);
    }
    i_start += 4;
  }
  // decode remaining height stripes
  height = block->size.y % 4;
  j_start = 0;
  for (uint16_t n2 = 0; n2 < num_h_stripe; n2++) {
    process_stripes_block(SigProp, block, i_start, j_start, width, height,
                          dum_stride, pLSB);
    j_start += 4;
  }
  width_last = block->size.x % 4;
  if (width_last) {
    process_stripes_block(SigProp, block, i_start, j_start, width_last, height,
                          dum_stride, pLSB);
  }
}

void ht_magref_decode(j2k_codeblock *block, uint8_t *HT_magref_segment,
                      uint32_t magref_length, const uint8_t &pLSB) {
  MR_dec MagRef(HT_magref_segment, magref_length);
  const uint16_t blk_height = block->size.y;
  const uint16_t blk_width = block->size.x;
  const uint16_t num_v_stripe = block->size.y / 4;
  uint16_t i_start = 0;
  uint16_t height = 4;
  int32_t *sp;

  for (uint16_t n1 = 0; n1 < num_v_stripe; n1++) {
    for (uint16_t j = 0; j < blk_width; j++) {
      for (uint16_t i = i_start; i < i_start + height; i++) {
        sp = &block->sample_buf[j + i * block->size.x];
        if (block->get_state(Sigma, i, j) != 0) {
          block->modify_state(refinement_indicator, 1, i, j);
          sp[0] |= MagRef.importMagRefBit() << pLSB;
        }
      }
    }
    i_start += 4;
  }
  height = blk_height % 4;
  for (uint16_t j = 0; j < blk_width; j++) {
    for (uint16_t i = i_start; i < i_start + height; i++) {
      sp = &block->sample_buf[j + i * block->size.x];
      if (block->get_state(Sigma, i, j) != 0) {
        block->modify_state(refinement_indicator, 1, i, j);
        sp[0] |= MagRef.importMagRefBit() << pLSB;
      }
    }
  }
}

int32_t htj2k_encode(j2k_codeblock *const block,
                     const uint8_t ROIshift) noexcept {
  /**
   write code here
   * */
}

int main() { return 0; }