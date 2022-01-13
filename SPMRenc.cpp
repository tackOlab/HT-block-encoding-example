#include "SPMRenc.h"
#include "inputs.h"
#include "codingunits.h"

void j2k_codeblock::set_MagSgn_and_sigma(uint32_t &or_val) {
  const uint32_t height = this->size.y;
  const uint32_t width  = this->size.x;
  const uint32_t stride = this->band_stride;

  int32_t *pp = in_sigma;
  for (uint16_t i = 0; i < height; ++i) {
    sprec_t *const sp  = this->i_samples + i * stride;
    int32_t *const dp  = this->sample_buf.get() + i * width;
    size_t block_index = (i + 1) * (size.x + 2) + 1;
    for (uint16_t j = 0; j < width; ++j) {
      int32_t temp  = sp[j];
      uint32_t sign = static_cast<uint32_t>(temp) & 0x80000000;
      block_states[block_index] |= *pp++;
      if (temp) {
        or_val |= 1;
      }
      //   block_states[block_index] |= 1;
      //   // convert sample value to MagSgn
      //   temp = (temp < 0) ? -temp : temp;
      //   temp &= 0x7FFFFFFF;
      //   temp--;
      //   temp <<= 1;
      //   temp += sign >> 31;
      dp[j] = temp;
      // }

      // if (sign) {
      //   printf("-%d ", dp[j] & 0x7fffffff);
      // } else {
      //   printf(" %d ", dp[j] & 0x7fffffff);
      // }
      block_index++;
    }
    // printf("\n");
  }
}

auto process_stripes_block = [](SP_enc &SigProp, j2k_codeblock *block, const uint16_t i_start,
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
      sp          = &block->sample_buf[j + i * block->size.x];
      causal_cond = (((block->Cmodes & CAUSAL) == 0) || (i != i_start + height - 1));
      mbr         = 0;
      if (block->get_state(Sigma, i, j) == 0) {
        block->calc_mbr(mbr, i, j, mbr_info & 0x1EF, causal_cond);
      }
      mbr_info >>= 3;
      if (mbr != 0) {
        bit = *sp & (1 << pLSB);  // HT Cleanupの設定によってpLSBは変化する可能性がある（当分0）
        SigProp.emitSPBit(bit);
        block->modify_state(refinement_indicator, 1, i, j);
        block->modify_state(refinement_value, bit, i, j);
      }
      block->modify_state(scan, 1, i, j);
    }
  }
  for (uint16_t j = j_start; j < j_start + width; j++) {
    for (uint16_t i = i_start; i < i_start + height; i++) {
      sp = &block->sample_buf[j + i * block->size.x];
      // encode sign
      // if ((((*sp & 0x7fffffff) & (1 << pLSB)) != 0) & (block->get_state(Refinement_indicator, i, j))) {
      if (block->get_state(Refinement_value, i, j)) {
        bit = (*sp < 0);  // sp < 0 なら bit = 1 となる
        SigProp.emitSPBit(bit);
      }
    }
  }
};

void ht_sigprop_encode(j2k_codeblock *block, SP_enc &SigProp, const uint8_t &pLSB) {
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
      process_stripes_block(SigProp, block, i_start, j_start, width, height, dum_stride, pLSB);
      j_start += 4;
    }
    width_last = block->size.x % 4;
    if (width_last) {
      process_stripes_block(SigProp, block, i_start, j_start, width_last, height, dum_stride, pLSB);
    }
    i_start += 4;
  }
  // encode remaining height stripes
  height  = block->size.y % 4;
  j_start = 0;
  for (uint16_t n2 = 0; n2 < num_h_stripe; n2++) {
    process_stripes_block(SigProp, block, i_start, j_start, width, height, dum_stride, pLSB);
    j_start += 4;
  }
  width_last = block->size.x % 4;
  if (width_last) {
    process_stripes_block(SigProp, block, i_start, j_start, width_last, height, dum_stride, pLSB);
  }
}

void ht_magref_encode(j2k_codeblock *block, MR_enc &MagRef, const uint8_t &pLSB) {
  const uint16_t blk_height   = block->size.y;
  const uint16_t blk_width    = block->size.x;
  const uint16_t num_v_stripe = block->size.y / 4;
  uint16_t i_start            = 0;
  uint16_t height             = 4;
  int32_t *sp;
  uint8_t bit;

  for (uint16_t n1 = 0; n1 < num_v_stripe; n1++) {
    for (uint16_t j = 0; j < blk_width; j++) {
      for (uint16_t i = i_start; i < i_start + height; i++) {
        sp = &block->sample_buf[j + i * block->size.x];
        if (block->get_state(Sigma, i, j) != 0) {
          bit = sp[0] & (1 << pLSB);  // HT Cleanupの設定によってpLSBは変化する可能性がある（当分0）
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
      sp = &block->sample_buf[j + i * block->size.x];
      if (block->get_state(Sigma, i, j) != 0) {
        bit = sp[0] & (1 << pLSB);  // HT Cleanupの設定によってpLSBは変化する可能性がある（当分0）
        MagRef.emitMRBit(bit);
        block->modify_state(refinement_indicator, 1, i, j);
      }
    }
  }
}

int32_t htj2k_encode(j2k_codeblock *const block, const uint8_t ROIshift) noexcept {
  uint32_t or_val = 0;
  // DWT係数をメンバ変数 sample_buf にコピー
  block->set_MagSgn_and_sigma(or_val);
  // SigProp + MagRef <= 2047 byte という決まりがある
  uint8_t Dref0[2047] = {0};
  uint8_t Dref1[2047] = {0};
  // エンコード状態のためのクラスを作成(initXXPacker)
  SP_enc SigProp(Dref0);
  MR_enc MagRef(Dref1);
  // SigProp encoding
  ht_sigprop_encode(block, SigProp, 0);
  // MagRef encoding
  ht_magref_encode(block, MagRef, 0);
  // 終端処理
  // MRのデータは左右反転してSPの後ろにくっつく
  termSPandMR(SigProp, MagRef);

  // 結果の表示
  printf("SP length = %d\n", SigProp.getPos());
  printf("MR length = %d\n", MagRef.getPos());
  for (int i = 0; i < SigProp.getPos() + MagRef.getPos(); ++i) {
    printf("%02X ", Dref0[i]);
  }
  printf("\n");
  return 0;
}

int main() {
  // コードブロックのパラメータ設定（今回だけ）
  uint32_t idx            = 0;
  uint8_t orientation     = 0;
  uint8_t M_b             = 8;
  uint8_t R_b             = 8;
  uint8_t transformation  = 0;
  float stepsize          = 1.0;
  uint32_t band_stride    = 64;
  uint32_t offset         = 0;
  uint16_t numlayers      = 1;
  uint8_t codeblock_style = 64;
  element_siz p0(0, 0);
  element_siz p1(64, 64);
  element_siz s(64, 64);

  // 量子化されたDWT係数を格納するバッファ領域の作成
  int32_t ibuf[4096];
  for (int i = 0; i < 4096; ++i) {
    ibuf[i] = (in_sign[i] << 31) | bitplane[i];
  }
  // コードブロッククラスのインスタンスの作成
  j2k_codeblock block(idx, orientation, M_b, R_b, transformation, stepsize, band_stride, ibuf, offset,
                      numlayers, codeblock_style, p0, p1, s);
  // エンコード処理
  htj2k_encode(&block, 0);

  return 0;
}