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

#include <memory>

#include "theta_sketch.hpp"

namespace datasketches {

/*
 * Updatable Theta sketch that operates in a contiguous region of memory of a fixed size.
 */
template<typename Allocator = std::allocator<uint64_t>>
class inplace_update_theta_sketch_alloc {
public:
  using Entry = uint64_t;
  using ExtractKey = trivial_extract_key;
  using Base = theta_update_sketch_base<Entry, ExtractKey, Allocator>;
  using iterator = theta_iterator<Entry, ExtractKey>;
  using const_iterator = theta_const_iterator<Entry, ExtractKey>;

  static void initialize(char* ptr, uint8_t lg_k, float p = 1, uint64_t seed = DEFAULT_SEED);

  explicit inplace_update_theta_sketch_alloc(char* ptr);

  static size_t size_bytes(uint8_t lg_k);
  static size_t size_u64(uint8_t lg_k);

  void update(uint64_t value);
  void update(const std::string& value);
  void update(const void* data, size_t length);

  void merge(const char* ptr);
  void merge_compact(const char* ptr);

  compact_theta_sketch_alloc<Allocator> compact(bool ordered = true, const Allocator& allocator = Allocator()) const;

  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;

private:
  struct inplace_update_theta_sketch_state {
    bool is_empty;
    uint8_t lg_k;
    uint32_t num_entries;
    uint64_t theta;
    uint64_t seed;
    uint64_t first_entry;
  };

  inplace_update_theta_sketch_state* state;

  // offsets are in sizeof(type)
  static const size_t COMPACT_SKETCH_PRE_LONGS_BYTE = 0;
  static const size_t COMPACT_SKETCH_FLAGS_BYTE = 5;
  static const size_t COMPACT_SKETCH_SEED_HASH_U16 = 3;
  static const size_t COMPACT_SKETCH_NUM_ENTRIES_U32 = 2;
  static const size_t COMPACT_SKETCH_SINGLE_ENTRY_U64 = 1;
  static const size_t COMPACT_SKETCH_ENTRIES_EXACT_U64 = 2;
  static const size_t COMPACT_SKETCH_THETA_U64 = 2;
  static const size_t COMPACT_SKETCH_ENTRIES_ESTIMATION_U64 = 3;

  static const uint8_t COMPACT_SKETCH_IS_EMPTY_FLAG = 2;

  void insert_or_ignore(uint64_t hash);
  void rebuild();
};

// alias with standard allocator for convenience
using inplace_update_theta_sketch = inplace_update_theta_sketch_alloc<std::allocator<uint64_t>>;

} /* namespace datasketches */

#include "inplace_update_theta_sketch_impl.hpp"

#endif