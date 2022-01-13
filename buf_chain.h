#pragma once
#include <cassert>
#include <cstdint>
#include <memory>
#include <vector>
/********************************************************************************
 * buf_chain
 *******************************************************************************/
class buf_chain {
 private:
  uint32_t total_length;
  uint32_t num_nodes;
  std::vector<uint8_t *> node_buf;
  std::vector<uint32_t> node_length;
  uint32_t node_pos;
  uint32_t pos;
  uint8_t *current_buf;
  uint32_t current_length;
  uint8_t tmp_byte;
  uint8_t last_byte;
  uint8_t bits;

 public:
  buf_chain() {
    total_length = 0;
    num_nodes = 0;
    node_buf = {};
    node_length = {};
    node_pos = 0;
    pos = 0;
    current_buf = nullptr;
    current_length = 0;
    tmp_byte = 0;
    last_byte = 0;
    bits = 0;
  }
  explicit buf_chain(uint32_t num) {
    total_length = 0;
    num_nodes = num;
    for (uint32_t i = 0; i < num; ++i) {
      node_buf.push_back(nullptr);
      node_length.push_back(0);
    }
    node_pos = 0;
    pos = 0;
    current_buf = nullptr;
    current_length = 0;
    tmp_byte = 0;
    last_byte = 0;
    bits = 0;
  }
  buf_chain &operator=(const buf_chain &bc) {
    if (this != &bc) {
      this->total_length = bc.total_length;
      this->num_nodes = bc.num_nodes;
      this->node_pos = bc.node_pos;
      this->pos = bc.pos;
      this->current_length = bc.current_length;
      this->tmp_byte = bc.tmp_byte;
      this->last_byte = bc.last_byte;
      this->bits = bc.bits;
      this->node_buf.reserve(bc.node_buf.size());
      for (size_t i = 0; i < bc.node_buf.size(); ++i) {
        this->node_buf.push_back(bc.node_buf[i]);
        this->node_length.push_back(bc.node_length[i]);
      }
    }
    return *this;
  }
  void add_buf_node(uint8_t *buf, uint32_t len) {
    total_length += len;
    node_buf.push_back(buf);
    node_length.push_back(len);
    num_nodes++;
  }
  void set_buf_node(uint32_t index, uint8_t *buf, uint32_t len) {
    total_length += len;
    node_buf[index] = buf;
    node_length[index] = len;
  }
  void activate() {
    pos = 0;
    node_pos = 0;
    current_buf = node_buf[0];
    current_length = node_length[0];
  }
  void activate(int n) {
    assert(n < this->node_buf.size());
    pos = 0;
    node_pos = n;
    current_buf = node_buf[node_pos];
    current_length = node_length[node_pos];
  }
  void flush_bits() { bits = 0; }
  void check_last_FF() {
    if (tmp_byte == 0xFF) {
      this->get_bit();
    }
  }
  uint32_t get_total_length() const { return total_length; }
  uint8_t get_specific_byte(uint32_t bufpos) { return *(current_buf + bufpos); }
  uint8_t get_byte() {
    if (pos > current_length - 1) {
      node_pos++;
      assert(node_pos <= num_nodes);
      current_buf = node_buf[node_pos];
      current_length = node_length[node_pos];
      pos = 0;
    }
    return *(current_buf + pos++);
  }
  uint8_t *get_current_address() {
    if (pos > current_length - 1) {
      node_pos++;
      assert(node_pos <= num_nodes);
      current_buf = node_buf[node_pos];
      current_length = node_length[node_pos];
      pos = 0;
    }
    return (current_buf + pos++);
  }
  void copy_N_bytes(uint8_t *&buf, uint32_t N) {
    // The first input argument is a reference of pointer.
    // The reason of this is the address of 'buf' shall be increased by N.
    assert((pos + N) <= current_length);
    memcpy(buf, current_buf + pos, N);
    pos += N;
    buf += N;
  }
  uint16_t get_word() {
    uint16_t word = get_byte();
    word <<= 8;
    word += get_byte();
    return word;
  }
  uint8_t get_bit() {
    if (bits == 0) {
      tmp_byte = get_byte();
      if (last_byte == 255) {
        bits = 7;
      } else {
        bits = 8;
      }
      last_byte = tmp_byte;
    }
    bits--;
    return (tmp_byte >> bits) & 1;
  }
  uint32_t get_N_bits(uint8_t N) {
    uint32_t cwd = 0;
    uint8_t bit;
    for (int i = 0; i < N; i++) {
      bit = get_bit();
      cwd <<= 1;
      cwd += static_cast<uint32_t>(bit);
    }
    return cwd;
  }
};
