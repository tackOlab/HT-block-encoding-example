#include "codingunits.h"
void j2k_codeblock::calc_mbr(uint8_t &mbr, uint16_t i, uint16_t j, uint32_t mbr_info,
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

  mbr |= get_state(Refinement_value, i - 1, j - 1) * get_state(Scan, i - 1, j - 1);
  mbr |= get_state(Refinement_value, i - 1, j) * get_state(Scan, i - 1, j);
  mbr |= get_state(Refinement_value, i - 1, j + 1) * get_state(Scan, i - 1, j + 1);
  mbr |= get_state(Refinement_value, i, j - 1) * get_state(Scan, i, j - 1);
  mbr |= get_state(Refinement_value, i, j + 1) * get_state(Scan, i, j + 1);
  mbr |= get_state(Refinement_value, i + 1, j - 1) * get_state(Scan, i + 1, j - 1) * causal_cond;
  mbr |= get_state(Refinement_value, i + 1, j) * get_state(Scan, i + 1, j) * causal_cond;
  mbr |= get_state(Refinement_value, i + 1, j + 1) * get_state(Scan, i + 1, j + 1) * causal_cond;
}

/********************************************************************************
 * j2k_codeblock
 *******************************************************************************/

j2k_codeblock::j2k_codeblock(const uint32_t &idx, uint8_t orientation, uint8_t M_b, uint8_t R_b,
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
      already_included(false),
      refsegment(false) {
  memset(sample_buf.get(), 0, sizeof(int32_t) * size.x * size.y);
  memset(block_states.get(), 0, (size.x + 2) * (size.y + 2));
  this->layer_start  = std::make_unique<uint8_t[]>(num_layers);
  this->layer_passes = std::make_unique<uint8_t[]>(num_layers);
  this->pass_length.reserve(109);
  this->pass_length = std::vector<uint32_t>(num_layers, 0);  // critical section
}

uint8_t j2k_codeblock::get_Mb() const { return this->M_b; }

uint8_t *j2k_codeblock::get_compressed_data() { return this->compressed_data.get(); }

void j2k_codeblock::set_compressed_data(uint8_t *const buf, const uint16_t bufsize, const uint16_t Lref) {
  if (this->compressed_data != nullptr) {
    if (!refsegment) {
      printf(
          "ERROR: illegal attempt to allocate codeblock's compressed data but the data is not "
          "null.\n");
      throw std::exception();
    } else {
      // if we are here, this function has been called to copy Dref[]
      memcpy(this->current_address + this->pass_length[0], buf, bufsize);
      return;
    }
  }
  this->compressed_data = std::make_unique<uint8_t[]>(bufsize + Lref * (refsegment));
  memcpy(this->compressed_data.get(), buf, bufsize);
  this->current_address = this->compressed_data.get();
}

void j2k_codeblock::create_compressed_buffer(buf_chain *tile_buf, uint16_t buf_limit,
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
        memcpy(this->compressed_data.get(), old_buf, sizeof(uint8_t) * (buf_limit));
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
