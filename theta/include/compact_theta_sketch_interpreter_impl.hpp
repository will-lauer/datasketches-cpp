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

#ifndef COMPACT_THETA_SKETCH_INTERPRETER_IMPL_HPP_
#define COMPACT_THETA_SKETCH_INTERPRETER_IMPL_HPP_

#include <iostream>
#include <iomanip>

namespace datasketches {

template<bool dummy>
auto compact_theta_sketch_interpreter<dummy>::interpret(const char* ptr, size_t size, uint64_t seed) -> compact_theta_sketch_data {
  if (size < 8) throw std::invalid_argument("at least 8 bytes expected, actual " + std::to_string(size));
  checker<true>::check_serial_version(ptr[COMPACT_SKETCH_SERIAL_VERSION_BYTE], COMPACT_SKETCH_SERIAL_VERSION);
  checker<true>::check_sketch_type(ptr[COMPACT_SKETCH_TYPE_BYTE], COMPACT_SKETCH_TYPE);
  uint64_t theta = theta_constants::MAX_THETA;
  if (ptr[COMPACT_SKETCH_FLAGS_BYTE] & (1 << COMPACT_SKETCH_IS_EMPTY_FLAG)) return {true, true, 0, theta, nullptr};
  const uint16_t seed_hash = reinterpret_cast<const uint16_t*>(ptr)[COMPACT_SKETCH_SEED_HASH_U16];
  checker<true>::check_seed_hash(seed_hash, compute_seed_hash(seed));
  const bool has_theta = ptr[COMPACT_SKETCH_PRE_LONGS_BYTE] > 2;
  if (has_theta) {
    if (size < 16) throw std::invalid_argument("at least 16 bytes expected, actual " + std::to_string(size));
    theta = reinterpret_cast<const uint64_t*>(ptr)[COMPACT_SKETCH_THETA_U64];
  }
  if (ptr[COMPACT_SKETCH_PRE_LONGS_BYTE] == 1) {
    return {false, true, 1, theta, reinterpret_cast<const uint64_t*>(ptr) + COMPACT_SKETCH_SINGLE_ENTRY_U64};
  }
  const uint32_t num_entries = reinterpret_cast<const uint32_t*>(ptr)[COMPACT_SKETCH_NUM_ENTRIES_U32];
  const size_t entries_start_u64 = has_theta ? COMPACT_SKETCH_ENTRIES_ESTIMATION_U64 : COMPACT_SKETCH_ENTRIES_EXACT_U64;
  const uint64_t* entries = reinterpret_cast<const uint64_t*>(ptr) + entries_start_u64;
  const size_t expected_size_bytes = (entries_start_u64 + num_entries) * sizeof(uint64_t);
  if (size < expected_size_bytes)
    throw std::invalid_argument(std::to_string(expected_size_bytes)
        + " bytes expected, actual " + std::to_string(size) + ", sketch dump: " + hex_dump(ptr, size));
  const bool is_ordered = ptr[COMPACT_SKETCH_FLAGS_BYTE] & (1 << COMPACT_SKETCH_IS_ORDERED_FLAG);
  return {false, is_ordered, num_entries, theta, entries};
}

template<bool dummy>
std::string compact_theta_sketch_interpreter<dummy>::hex_dump(const char* ptr, size_t size) {
  std::stringstream s;
  s << std::hex << std::setfill('0') << std::uppercase;
  for (size_t i = 0; i < size; ++i) s << std::setw(2) << static_cast<int>(ptr[i]);
  return s.str();
}

} /* namespace datasketches */

#endif
