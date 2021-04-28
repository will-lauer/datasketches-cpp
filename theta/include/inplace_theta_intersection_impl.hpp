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

#ifndef INPLACE_THETA_INTERSECTION_IMPL_HPP_
#define INPLACE_THETA_INTERSECTION_IMPL_HPP_

#include "common_defs.hpp"
#include "theta_helpers.hpp"
#include "compact_theta_sketch_interpreter.hpp"

namespace datasketches {

template<typename B, typename A>
void inplace_theta_intersection_alloc<B, A>::initialize(B& buffer, uint64_t seed) {
  buffer.resize(sizeof(uint64_t) * 2); // store seed here until the first update and invalid marker of the same size for alignment
  auto data = reinterpret_cast<uint64_t*>(buffer.data());
  data[0] = INVALID_MARKER; // invalid state marker, impossible for a valid sketch
  data[1] = seed; // store seed here until table is needed
}

template<typename B, typename A>
inplace_theta_intersection_alloc<B, A>::inplace_theta_intersection_alloc(B& buffer):
buffer_(buffer) {}

template<typename B, typename A>
size_t inplace_theta_intersection_alloc<B, A>::max_size_bytes(uint8_t lg_k) {
  return Sketch::max_size_bytes(lg_k);
}

template<typename B, typename A>
void inplace_theta_intersection_alloc<B, A>::intersection(const char* ptr, size_t size) {
  if (size < sizeof(uint64_t)) throw std::invalid_argument("at least " + std::to_string(sizeof(uint64_t))
      + " bytes expected, actual " + std::to_string(size));
  if (reinterpret_cast<const uint64_t*>(ptr)[0] == INVALID_MARKER) throw std::invalid_argument("invalid input");
  auto other_state = reinterpret_cast<const typename Sketch::inplace_update_theta_sketch_state*>(ptr);
  const size_t expected_size_bytes = Sketch::header_size_bytes() + Sketch::table_size_bytes(other_state->lg_cur_size);
  if (size < expected_size_bytes) throw std::invalid_argument(std::to_string(expected_size_bytes)
        + " bytes expected, actual " + std::to_string(size));
  auto data = reinterpret_cast<uint64_t*>(buffer_.data());
  if (data[0] == INVALID_MARKER) { // first update, create sketch and insert entries from the input sketch
    uint64_t seed = data[1];
    if (seed != other_state->seed) throw std::invalid_argument("seed mismatch");
    Sketch::initialize(buffer_, other_state->lg_cur_size, other_state->lg_cur_size, Sketch::resize_factor::X1, other_state->theta, seed, other_state->is_empty);
    Sketch sketch(buffer_);
    const Entry* other_entries = &other_state->first_entry;
    const uint32_t other_size = 1 << other_state->lg_cur_size;
    for (uint32_t i = 0; i < other_size; ++i) if (other_entries[i]) sketch.insert_or_ignore(other_entries[i]);
  } else { // intersection
    Sketch sketch(buffer_);
    if (sketch.is_empty()) return;
    if (sketch.seed() != other_state->seed) throw std::invalid_argument("seed mismatch");
    sketch.is_empty(sketch.is_empty() | other_state->is_empty);
    sketch.theta(std::min(sketch.theta(), other_state->theta));
    if (sketch.num_entries() == 0) return;
    if (other_state->num_entries == 0) {
      Sketch::initialize(buffer_, 0, 0, Sketch::resize_factor::X1, sketch.theta(), sketch.seed(), sketch.is_empty());
      return;
    }

    const uint32_t max_matches = std::min(sketch.num_entries(), other_state->num_entries);
    std::vector<uint64_t, A> matched_entries;
    matched_entries.reserve(max_matches);
    uint32_t count = 0;
    const Entry* other_entries = &other_state->first_entry;
    const uint32_t other_size = 1 << other_state->lg_cur_size;
    for (uint32_t i = 0; i < other_size; ++i) {
      const uint64_t entry = other_entries[i];
      if (entry != 0) {
        if(entry < sketch.theta()) {
          if (sketch.contains(entry)) {
            if (matched_entries.size() == max_matches) throw std::invalid_argument("max matches exceeded, possibly corrupted input sketch");
            matched_entries.push_back(entry);
          }
        }
        ++count;
      }
    }
    if (count != other_state->num_entries) {
      throw std::invalid_argument(std::to_string(other_state->num_entries) + " keys expected, actual "
          + std::to_string(count) + ", possibly corrupted input sketch");
    }
    if (matched_entries.size() == 0) {
      Sketch::initialize(buffer_, 0, 0, Sketch::resize_factor::X1, sketch.theta(), sketch.seed(), sketch.is_empty());
      if (sketch.theta() == theta_constants::MAX_THETA) sketch.is_empty(true);
    } else {
      const uint8_t lg_size = lg_size_from_count(matched_entries.size(), Base::REBUILD_THRESHOLD);
      Sketch::initialize(buffer_, lg_size, lg_size, Sketch::resize_factor::X1, sketch.theta(), sketch.seed(), sketch.is_empty());
      for (auto entry: matched_entries) sketch.insert_or_ignore(entry);
    }
  }
}

template<typename B, typename A>
void inplace_theta_intersection_alloc<B, A>::intersection_compact(const char* ptr, size_t size) {
  auto data = reinterpret_cast<uint64_t*>(buffer_.data());
  if (data[0] == INVALID_MARKER) { // first update, create sketch and insert entries from the input sketch
    uint64_t seed = data[1];
    auto input_sketch_data = compact_theta_sketch_interpreter<true>::interpret(ptr, size, seed);
    const uint8_t lg_size = lg_size_from_count(input_sketch_data.num_entries, Base::REBUILD_THRESHOLD);
    Sketch::initialize(buffer_, lg_size, lg_size, Sketch::resize_factor::X1, input_sketch_data.theta, seed, input_sketch_data.is_empty);
    Sketch sketch(buffer_);
    for (size_t i = 0; i < input_sketch_data.num_entries; ++i) sketch.insert_or_ignore(input_sketch_data.entries[i]);
  } else { // intersection
    Sketch sketch(buffer_);
    if (sketch.is_empty()) return;
    auto input_sketch_data = compact_theta_sketch_interpreter<true>::interpret(ptr, size, sketch.seed());
    sketch.is_empty(sketch.is_empty() | input_sketch_data.is_empty);
    sketch.theta(std::min(sketch.theta(), input_sketch_data.theta));
    if (sketch.num_entries() == 0) return;
    if (input_sketch_data.num_entries == 0) {
      Sketch::initialize(buffer_, 0, 0, Sketch::resize_factor::X1, sketch.theta(), sketch.seed(), sketch.is_empty());
      return;
    }

    const uint32_t max_matches = std::min(sketch.num_entries(), input_sketch_data.num_entries);
    std::vector<uint64_t, A> matched_entries;
    matched_entries.reserve(max_matches);
    uint32_t count = 0;
    for (uint32_t i = 0; i < input_sketch_data.num_entries; ++i) {
      const uint64_t entry = input_sketch_data.entries[i];
      if (entry < sketch.theta()) {
        if (sketch.contains(entry)) {
          if (matched_entries.size() == max_matches) throw std::invalid_argument("max matches exceeded, possibly corrupted input sketch");
          matched_entries.push_back(entry);
        }
      } else if (input_sketch_data.is_ordered) {
        break; // early stop
      }
      ++count;
    }
    if (count > input_sketch_data.num_entries) {
      throw std::invalid_argument("more keys than expected, possibly corrupted input sketch");
    } else if (!input_sketch_data.is_ordered && count < input_sketch_data.num_entries) {
      throw std::invalid_argument("fewer keys than expected, possibly corrupted input sketch");
    }
    if (matched_entries.size() == 0) {
      Sketch::initialize(buffer_, 0, 0, Sketch::resize_factor::X1, sketch.theta(), sketch.seed(), sketch.is_empty());
      if (sketch.theta() == theta_constants::MAX_THETA) sketch.is_empty(true);
    } else {
      const uint8_t lg_size = lg_size_from_count(matched_entries.size(), Base::REBUILD_THRESHOLD);
      Sketch::initialize(buffer_, lg_size, lg_size, Sketch::resize_factor::X1, sketch.theta(), sketch.seed(), sketch.is_empty());
      for (auto entry: matched_entries) sketch.insert_or_ignore(entry);
    }
  }
}

template<typename B, typename A>
compact_theta_sketch_alloc<A> inplace_theta_intersection_alloc<B, A>::get_result(bool ordered, const A& allocator) const {
  auto data = reinterpret_cast<uint64_t*>(buffer_.data());
  if (data[0] == INVALID_MARKER) throw std::invalid_argument("calling get_result() before the first intersection is undefined");
  return Sketch(buffer_).compact(ordered, allocator);
}

} /* namespace datasketches */

#endif
