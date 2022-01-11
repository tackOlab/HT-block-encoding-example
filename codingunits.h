#include <cstdio>
#include <cstring>
#include <functional>

#include "buf_chain.h"

typedef int32_t sprec_t;

class element_siz {
 public:
  uint32_t x;
  uint32_t y;
  element_siz() : x(0), y(0) {}
  element_siz(uint32_t x0, uint32_t y0) {
    x = x0;
    y = y0;
  }
};

#define CAUSAL 0x008

#define SHIFT_SIGMA 0   // J2K and HTJ2K
#define SHIFT_SIGMA_ 1  // J2K only
#define SHIFT_PI_ 2     // J2K and HTJ2K; used as refinement indicator for HTJ2K
#define SHIFT_REF 3     // HTJ2K only
#define SHIFT_SCAN 4    // HTJ2K only
#define SHIFT_P 3       // J2K only

// getters
inline uint8_t Sigma(uint8_t &data) { return (data >> SHIFT_SIGMA) & 1; }
inline uint8_t Sigma_(uint8_t &data) { return (data >> SHIFT_SIGMA_) & 1; }
inline uint8_t Pi_(uint8_t &data) { return (data >> SHIFT_PI_) & 1; }
inline uint8_t Scan(uint8_t &data) { return (data >> SHIFT_SCAN) & 1; }
inline uint8_t Refinement_value(uint8_t &data) {
  return (data >> SHIFT_REF) & 1;
}
inline uint8_t Refinement_indicator(uint8_t &data) {
  return (data >> SHIFT_PI_) & 1;
}
inline uint8_t Decoded_bitplane_index(uint8_t &data) {
  return (data >> SHIFT_P);
}

// setters
inline void sigma(uint8_t &data, const uint8_t &val) { data |= val; }
inline void sigma_(uint8_t &data, const uint8_t &val) {
  data |= val << SHIFT_SIGMA_;
}
inline void pi_(uint8_t &data, const uint8_t &val) {
  if (val) {
    data |= 1 << SHIFT_PI_;
  } else {
    data &= ~(1 << SHIFT_PI_);
  }
}
inline void scan(uint8_t &data, const uint8_t &val) {
  data |= val << SHIFT_SCAN;
}
inline void refinement_value(uint8_t &data, const uint8_t &val) {
  data |= val << SHIFT_REF;
}
inline void refinement_indicator(uint8_t &data, const uint8_t &val) {
  if (val) {
    data |= 1 << SHIFT_PI_;
  } else {
    data &= ~(1 << SHIFT_PI_);
  }
}
inline void decoded_bitplane_index(uint8_t &data, const uint8_t &val) {
  data &= 0x07;
  data |= val << SHIFT_P;
}

/********************************************************************************
 * j2k_region
 *******************************************************************************/
class j2k_region {
 public:
  // top-left coordinate (inclusive) of a region in the reference grid
  element_siz pos0;
  // bottom-right coordinate (exclusive) of a region in the reference grid
  element_siz pos1;
  // return top-left coordinate (inclusive)
  element_siz get_pos0() const { return pos0; }
  // return bottom-right coordinate (exclusive)
  element_siz get_pos1() const { return pos1; }
  // get size of a region
  void get_size(element_siz &out) const {
    out.x = pos1.x - pos0.x;
    out.y = pos1.y - pos0.y;
  }
  // set top-left coordinate (inclusive)
  void set_pos0(element_siz in) { pos0 = in; }
  // set bottom-right coordinate (exclusive)
  void set_pos1(element_siz in) { pos1 = in; }
  j2k_region() = default;
  j2k_region(element_siz p0, element_siz p1) : pos0(p0), pos1(p1) {}
};

/********************************************************************************
 * j2k_codeblock
 *******************************************************************************/
class j2k_codeblock : public j2k_region {
 public:
  const element_siz size;

 private:
  const uint32_t index;
  const uint8_t band;
  const uint8_t M_b;
  std::unique_ptr<uint8_t[]> compressed_data;
  uint8_t *current_address;

 public:
  std::unique_ptr<uint8_t[]> block_states;
  const uint8_t R_b;
  const uint8_t transformation;
  const float stepsize;
  const uint32_t band_stride;
  const uint16_t num_layers;
  std::unique_ptr<int32_t[]> sample_buf;
  sprec_t *const i_samples;
  uint32_t length;
  uint16_t Cmodes;
  uint8_t num_passes;
  uint8_t num_ZBP;
  uint8_t fast_skip_passes;
  uint32_t Lblock;
  // length of a coding pass in byte
  std::vector<uint32_t> pass_length;
  // index of the coding-pass from which layer starts
  std::unique_ptr<uint8_t[]> layer_start;
  // number of coding-passes included in a layer
  std::unique_ptr<uint8_t[]> layer_passes;
  bool already_included;

  j2k_codeblock(const uint32_t &idx, uint8_t orientation, uint8_t M_b,
                uint8_t R_b, uint8_t transformation, float stepsize,
                uint32_t band_stride, sprec_t *ibuf, uint32_t offset,
                const uint16_t &numlayers, const uint8_t &codeblock_style,
                const element_siz &p0, const element_siz &p1,
                const element_siz &s);
  void modify_state(const std::function<void(uint8_t &, uint8_t)> &callback,
                    uint8_t val, int16_t j1, int16_t j2) {
    callback(block_states[(j1 + 1) * (size.x + 2) + (j2 + 1)], val);
  }
  uint8_t get_state(const std::function<uint8_t(uint8_t &)> &callback,
                    int16_t j1, int16_t j2) const {
    return callback(block_states[(j1 + 1) * (size.x + 2) + (j2 + 1)]);
  }
  // DEBUG FUNCTION, SOON BE DELETED
  uint8_t get_orientation() const { return band; }
  uint8_t get_context_label_sig(const uint16_t &j1, const uint16_t &j2) const;
  uint8_t get_signLUT_index(const uint16_t &j1, const uint16_t &j2) const;
  uint8_t get_Mb() const;
  uint8_t *get_compressed_data();
  void set_compressed_data(uint8_t *buf, uint16_t size);
  void create_compressed_buffer(buf_chain *tile_buf, uint16_t buf_limit,
                                const uint16_t &layer);
  void update_sample(const uint8_t &symbol, const uint8_t &p,
                     const uint16_t &j1, const uint16_t &j2) const;
  void update_sign(const int8_t &val, const uint16_t &j1,
                   const uint16_t &j2) const;
  uint8_t get_sign(const uint16_t &j1, const uint16_t &j2) const;
  void set_MagSgn_and_sigma(uint32_t &or_val);
  void calc_mbr(uint8_t &mbr, uint16_t i, uint16_t j, uint32_t mbr_info,
                uint8_t causal_cond) const {
    // mbr |= (mbr_info != 0);
    mbr |= get_state(Sigma, i - 1, j - 1);
    mbr |= get_state(Sigma, i - 1, j);
    mbr |= get_state(Sigma, i - 1, j + 1);
    mbr |= get_state(Sigma, i, j - 1);
    mbr |= get_state(Sigma, i, j + 1);
    mbr |= get_state(Sigma, i + 1, j - 1) * causal_cond;
    mbr |= get_state(Sigma, i + 1, j) * causal_cond;
    mbr |= get_state(Sigma, i + 1, j + 1) * causal_cond;

    mbr |= get_state(Refinement_value, i - 1, j - 1) *
           get_state(Scan, i - 1, j - 1);
    mbr |= get_state(Refinement_value, i - 1, j) * get_state(Scan, i - 1, j);
    mbr |= get_state(Refinement_value, i - 1, j + 1) *
           get_state(Scan, i - 1, j + 1);
    mbr |= get_state(Refinement_value, i, j - 1) * get_state(Scan, i, j - 1);
    mbr |= get_state(Refinement_value, i, j + 1) * get_state(Scan, i, j + 1);
    mbr |= get_state(Refinement_value, i + 1, j - 1) *
           get_state(Scan, i + 1, j - 1) * causal_cond;
    mbr |= get_state(Refinement_value, i + 1, j) * get_state(Scan, i + 1, j) *
           causal_cond;
    mbr |= get_state(Refinement_value, i + 1, j + 1) *
           get_state(Scan, i + 1, j + 1) * causal_cond;
  }
};

/********************************************************************************
 * j2k_codeblock
 *******************************************************************************/

j2k_codeblock::j2k_codeblock(
    const uint32_t &idx, uint8_t orientation, uint8_t M_b, uint8_t R_b,
    uint8_t transformation, float stepsize, uint32_t band_stride, sprec_t *ibuf,
    uint32_t offset, const uint16_t &numlayers, const uint8_t &codeblock_style,
    const element_siz &p0, const element_siz &p1, const element_siz &s)
    : j2k_region(p0, p1),
      // public
      size(s),
      // private
      index(idx),
      band(orientation),
      M_b(M_b),
      compressed_data(nullptr),
      current_address(nullptr),
      block_states(std::make_unique<uint8_t[]>((size.x + 2) * (size.y + 2))),
      // public
      R_b(R_b),
      transformation(transformation),
      stepsize(stepsize),
      band_stride(band_stride),
      num_layers(numlayers),
      sample_buf(std::make_unique<int32_t[]>(size.x * size.y)),
      i_samples(ibuf + offset),
      length(0),
      Cmodes(codeblock_style),
      num_passes(0),
      num_ZBP(0),
      fast_skip_passes(0),
      Lblock(0),
      already_included(false) {
  memset(sample_buf.get(), 0, sizeof(int32_t) * size.x * size.y);
  memset(block_states.get(), 0, (size.x + 2) * (size.y + 2));
  this->layer_start = std::make_unique<uint8_t[]>(num_layers);
  this->layer_passes = std::make_unique<uint8_t[]>(num_layers);
  this->pass_length.reserve(109);
  this->pass_length = std::vector<uint32_t>(num_layers, 0);  // critical section
}

uint8_t j2k_codeblock::get_Mb() const { return this->M_b; }

uint8_t *j2k_codeblock::get_compressed_data() {
  return this->compressed_data.get();
}

void j2k_codeblock::set_compressed_data(uint8_t *buf, uint16_t bufsize) {
  if (this->compressed_data != nullptr) {
    printf(
        "ERROR: illegal attempt to allocate codeblock's compressed data but "
        "the data is not "
        "null.\n");
    throw std::exception();
  }
  this->compressed_data = std::make_unique<uint8_t[]>(bufsize);
  memcpy(this->compressed_data.get(), buf, bufsize);
  this->current_address = this->compressed_data.get();
}

void j2k_codeblock::create_compressed_buffer(buf_chain *tile_buf,
                                             uint16_t buf_limit,
                                             const uint16_t &layer) {
  uint32_t layer_length = 0;
  uint16_t l0, l1;
  if (this->layer_passes[layer] > 0) {
    l0 = this->layer_start[layer];
    l1 = l0 + this->layer_passes[layer];
    for (int i = l0; i < l1; i++) {
      layer_length += this->pass_length[i];
    }
    // allocate buffer one once for the first contributing layer
    if (this->compressed_data == nullptr) {
      this->compressed_data = std::make_unique<uint8_t[]>(buf_limit);
      this->current_address = this->compressed_data.get();
    }
    if (layer_length != 0) {
      while (this->length + layer_length > buf_limit) {
        // extend buffer size, if necessary
        uint8_t *old_buf = this->compressed_data.release();
        buf_limit += 8192;
        this->compressed_data = std::make_unique<uint8_t[]>(buf_limit);
        memcpy(this->compressed_data.get(), old_buf,
               sizeof(uint8_t) * (buf_limit));
        this->current_address = this->compressed_data.get() + (this->length);
        delete[] old_buf;
      }
      // we assume that the size of the compressed data is less than or equal to
      // that of buf_chain node.
      tile_buf->copy_N_bytes(this->current_address, layer_length);
      this->length += layer_length;
    }
  }
}
