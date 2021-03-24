/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

#ifndef INPLACE_UPDATE_THETA_SKETCH_HPP_
#define INPLACE_UPDATE_THETA_SKETCH_HPP_

#include <vector>

#include "theta_sketch.hpp"

namespace datasketches {

/*
 * Updatable Theta sketch that operates in a contiguous region of memory of a fixed size.
 */
template<typename Buffer = std::vector<char>, typename Allocator = std::allocator<uint64_t>>
class inplace_update_theta_sketch_alloc {
public:
  using Entry = uint64_t;
  using ExtractKey = trivial_extract_key;
  using Base = theta_update_sketch_base<Entry, ExtractKey, Allocator>;
  using iterator = theta_iterator<Entry, ExtractKey>;
  using const_iterator = theta_const_iterator<Entry, ExtractKey>;
  using resize_factor = theta_constants::resize_factor;

  class builder;

  explicit inplace_update_theta_sketch_alloc(Buffer& buffer);

  static size_t max_size_bytes(uint8_t lg_k);

  void update(uint64_t value);
  void update(const std::string& value);
  void update(const void* data, size_t length);

  void merge(const char* ptr, size_t size);
  void merge_compact(const char* ptr, size_t size);

  inplace_update_theta_sketch_alloc& trim();

  compact_theta_sketch_alloc<Allocator> compact(bool ordered = true, const Allocator& allocator = Allocator()) const;

  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;

private:
  struct inplace_update_theta_sketch_state {
    bool is_empty;
    uint8_t lg_cur_size;
    uint8_t lg_nom_size;
    resize_factor rf;
    uint32_t num_entries;
    uint64_t theta;
    uint64_t seed;
    uint64_t first_entry;
  };

  Buffer& buffer;

  // offsets are in sizeof(type)
  static const size_t COMPACT_SKETCH_PRE_LONGS_BYTE = 0;
  static const size_t COMPACT_SKETCH_SERIAL_VERSION_BYTE = 1;
  static const size_t COMPACT_SKETCH_TYPE_BYTE = 2;
  static const size_t COMPACT_SKETCH_FLAGS_BYTE = 5;
  static const size_t COMPACT_SKETCH_SEED_HASH_U16 = 3;
  static const size_t COMPACT_SKETCH_NUM_ENTRIES_U32 = 2;
  static const size_t COMPACT_SKETCH_SINGLE_ENTRY_U64 = 1;
  static const size_t COMPACT_SKETCH_ENTRIES_EXACT_U64 = 2;
  static const size_t COMPACT_SKETCH_THETA_U64 = 2;
  static const size_t COMPACT_SKETCH_ENTRIES_ESTIMATION_U64 = 3;

  static const uint8_t COMPACT_SKETCH_IS_EMPTY_FLAG = 2;

  static const uint8_t COMPACT_SKETCH_SERIAL_VERSION = 3;
  static const uint8_t COMPACT_SKETCH_TYPE = 3;

  void insert_or_ignore(uint64_t hash);
  void resize();
  void rebuild();

  static size_t header_size_bytes();
  static size_t table_size_bytes(uint8_t lg_k);

  // for builder
  static void initialize(Buffer& buffer, uint8_t lg_cur_size, uint8_t lg_nom_size, resize_factor rf, uint64_t theta, uint64_t seed);
};

template<typename Buffer, typename Allocator>
class inplace_update_theta_sketch_alloc<Buffer, Allocator>::builder: public theta_base_builder<builder, Allocator> {
public:
  builder();
  void initialize(Buffer& buffer) const;
};

// alias with default types for convenience
using inplace_update_theta_sketch = inplace_update_theta_sketch_alloc<>;

} /* namespace datasketches */

#include "inplace_update_theta_sketch_impl.hpp"

#endif
