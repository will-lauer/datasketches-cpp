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

#ifndef INPLACE_THETA_INTERSECTION_HPP_
#define INPLACE_THETA_INTERSECTION_HPP_

#include <vector>
#include <tuple>

#include "inplace_update_theta_sketch.hpp"

namespace datasketches {

/*
 * Theta intersection that operates in a contiguous region of memory represented as std::vector<char>.
 */
template<typename Buffer = std::vector<char>, typename Allocator = std::allocator<uint64_t>>
class inplace_theta_intersection_alloc {
public:
  using Entry = uint64_t;
  using ExtractKey = trivial_extract_key;
  using Base = theta_update_sketch_base<Entry, ExtractKey, Allocator>;
  using Sketch = inplace_update_theta_sketch_alloc<Buffer, Allocator>;

  static void initialize(Buffer& buffer, uint64_t seed = DEFAULT_SEED);
  explicit inplace_theta_intersection_alloc(Buffer& buffer);

  static size_t max_size_bytes(uint8_t lg_k);

  void intersection(const char* ptr, size_t size);
  void intersection_compact(const char* ptr, size_t size);

  compact_theta_sketch_alloc<Allocator> get_result(bool ordered = true, const Allocator& allocator = Allocator()) const;

private:
  Buffer& buffer_;

  static const uint64_t INVALID_MARKER = 0xaa55aa55aa55aa55;
};

// alias with default types for convenience
using inplace_theta_intersection = inplace_theta_intersection_alloc<>;

} /* namespace datasketches */

#include "inplace_theta_intersection_impl.hpp"

#endif
