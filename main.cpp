#include <cstdio>
#include <cstdint>
#include "codingunits.h"
#include "SPMRenc.h"
#include "inputs.h"

int32_t htj2k_cleanup_encode(j2k_codeblock *const block, const uint8_t ROIshift) noexcept;
void ht_sigprop_encode(j2k_codeblock *block, SP_enc &SigProp);
void ht_magref_encode(j2k_codeblock *block, MR_enc &MagRef);

int32_t htj2k_encode(j2k_codeblock *const block, const uint8_t ROIshift) noexcept {
  // true for generating HT Refinement segment
  block->refsegment = true;
  int32_t Lcup      = htj2k_cleanup_encode(block, ROIshift);
  printf("- Length of HT Cleanup segment = %d\n", Lcup);
  if (Lcup && block->refsegment) {
    // SigProp + MagRef <= 2047 byte という決まりがある
    uint8_t Dref[2047] = {0};
    // uint8_t Dref1[2047] = {0};
    //  エンコード状態のためのクラスを作成(initXXPacker)
    SP_enc SigProp(Dref);
    MR_enc MagRef(Dref);
    // SigProp encoding
    ht_sigprop_encode(block, SigProp);
    // MagRef encoding
    ht_magref_encode(block, MagRef);
    // 終端処理
    int32_t HTRefinement_segment = 0;
    if (MagRef.get_length()) {
      HTRefinement_segment = termSPandMR(SigProp, MagRef);
      block->num_passes += 2;
      block->layer_passes[0] += 2;
      block->pass_length.push_back(SigProp.get_length());
      block->pass_length.push_back(MagRef.get_length());
    } else {
      SigProp.termSP();
      HTRefinement_segment = SigProp.get_length();
      block->num_passes += 1;
      block->layer_passes[0] += 1;
      block->pass_length.push_back(SigProp.get_length());
    }
    if (HTRefinement_segment) {
      block->length += HTRefinement_segment;
      block->num_ZBP -= (block->refsegment);
      block->set_compressed_data(Dref, HTRefinement_segment);
    }
    // 結果の表示
    printf("- Length of HT Refinement segment = %d\n", HTRefinement_segment);
    printf("  - SP length = %d\n", SigProp.get_length());
    printf("  - MR length = %d\n", MagRef.get_length());
    // for (int i = 0; i < HTRefinement_segment; ++i) {
    //   printf("%02X ", Dref[i]);
    // }
    printf("\n");
  }
  printf("Number of coding passes = %d\n", block->num_passes);
  printf("Value of Z_blk = %d\n", block->num_ZBP);
  printf("Value of eps_b = %d\n", block->get_Mb() - 1);  // assuming guard bits = 1
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
  uint8_t codeblock_style = 0b01000000;
  element_siz p0(0, 0);
  element_siz p1(64, 64);
  element_siz s(64, 64);

  // コードブロッククラスのインスタンスの作成
  j2k_codeblock block(idx, orientation, M_b, R_b, transformation, stepsize, band_stride, quantized_samples,
                      offset, numlayers, codeblock_style, p0, p1, s);
  // エンコード処理
  htj2k_encode(&block, 0);

  return 0;
}