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

#ifndef INPLACE_UPDATE_THETA_SKETCH_IMPL_HPP_
#define INPLACE_UPDATE_THETA_SKETCH_IMPL_HPP_

#include "theta_helpers.hpp"

namespace datasketches {

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::initialize(B& buffer, uint8_t lg_cur_size, uint8_t lg_nom_size, resize_factor rf, uint64_t theta, uint64_t seed) {
  buffer.resize(header_size_bytes() + table_size_bytes(lg_cur_size));
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  state->is_empty = true;
  state->lg_cur_size = lg_cur_size;
  state->lg_nom_size = lg_nom_size;
  state->rf = rf;
  state->num_entries = 0;
  state->theta = theta;
  state->seed = seed;
  Entry* entries = &state->first_entry;
  std::fill(entries, entries + (1 << lg_cur_size), 0); // might be unnecessary
}

template<typename B, typename A>
inplace_update_theta_sketch_alloc<B, A>::inplace_update_theta_sketch_alloc(B& buffer):
buffer(buffer) {}

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::update(uint64_t value) {
  update(&value, sizeof(value));
}

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::update(const std::string& value) {
  if (value.empty()) return;
  update(value.c_str(), value.length());
}

template<typename B, typename A>
size_t inplace_update_theta_sketch_alloc<B, A>::max_size_bytes(uint8_t lg_k) {
  return header_size_bytes() + table_size_bytes(lg_k + 1);
}

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::update(const void* data, size_t length) {
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  state->is_empty = false;
  const uint64_t hash = compute_hash(data, length, state->seed);
  insert_or_ignore(hash);
}

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::insert_or_ignore(uint64_t hash) {
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  if (hash >= state->theta) return; // hash == 0 is reserved to mark empty slots in the table
  uint64_t* entries = &state->first_entry;
  auto result = Base::find(entries, state->lg_cur_size, hash);
  if (!result.second) {
    *result.first = hash;
    ++state->num_entries;
    if (state->num_entries > Base::get_capacity(state->lg_cur_size, state->lg_nom_size)) {
      if (state->lg_cur_size <= state->lg_nom_size) {
        resize();
      } else {
        rebuild();
      }
    }
  }
}

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::merge(const char* ptr) {
  auto other_state = reinterpret_cast<const inplace_update_theta_sketch_state*>(ptr);
  if (other_state->is_empty) return;
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  if (state->seed != other_state->seed) throw std::invalid_argument("seed mismatch");
  state->is_empty = false;
  if (state->theta > other_state->theta) {
    state->theta = other_state->theta;
    std::vector<uint64_t, A> new_entries;
    new_entries.reserve(state->num_entries);
    for (auto entry: *this) {
      if (entry < state->theta) new_entries.push_back(entry);
    }
    Entry* entries = &state->first_entry;
    std::fill(entries, entries + (1 << state->lg_cur_size), 0);
    for (auto entry: new_entries) {
      auto result = Base::find(entries, state->lg_cur_size, entry);
      if (!result.second) *result.first = entry;
    }
    state->num_entries = new_entries.size();
  }
  const Entry* other_entries = &other_state->first_entry;
  const size_t other_size = 1 << other_state->lg_cur_size;
  for (size_t i = 0; i < other_size; ++i) if (other_entries[i]) insert_or_ignore(other_entries[i]);
}

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::merge_compact(const char* ptr) {
  // TODO: check serial version and sketch type
  if (ptr[COMPACT_SKETCH_FLAGS_BYTE] & (1 << COMPACT_SKETCH_IS_EMPTY_FLAG)) return;
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  const uint16_t other_seed_hash = reinterpret_cast<const uint16_t*>(ptr)[COMPACT_SKETCH_SEED_HASH_U16];
  checker<true>::check_seed_hash(other_seed_hash, compute_seed_hash(state->seed));
  state->is_empty = false;
  const bool other_has_theta = ptr[COMPACT_SKETCH_PRE_LONGS_BYTE] > 2;
  if (other_has_theta) {
    const uint64_t other_theta = reinterpret_cast<const uint64_t*>(ptr)[COMPACT_SKETCH_THETA_U64];
    if (state->theta > other_theta) {
      state->theta = other_theta;
      // TODO: extract common rebuild code
      std::vector<uint64_t, A> new_entries;
      new_entries.reserve(state->num_entries);
      for (auto entry: *this) {
        if (entry < state->theta) new_entries.push_back(entry);
      }
      Entry* entries = &state->first_entry;
      std::fill(entries, entries + (1 << state->lg_cur_size), 0);
      for (auto entry: new_entries) {
        auto result = Base::find(entries, state->lg_cur_size, entry);
        if (!result.second) *result.first = entry;
      }
      state->num_entries = new_entries.size();
    }
  }
  if (ptr[COMPACT_SKETCH_PRE_LONGS_BYTE] == 1) {
    insert_or_ignore(reinterpret_cast<const uint64_t*>(ptr)[COMPACT_SKETCH_SINGLE_ENTRY_U64]);
  } else {
    const size_t num_entries = reinterpret_cast<const uint32_t*>(ptr)[COMPACT_SKETCH_NUM_ENTRIES_U32];
    const Entry* other_entries = reinterpret_cast<const uint64_t*>(ptr) +
        (other_has_theta ? COMPACT_SKETCH_ENTRIES_ESTIMATION_U64 : COMPACT_SKETCH_ENTRIES_EXACT_U64);
    for (size_t i = 0; i < num_entries; ++i) insert_or_ignore(other_entries[i]);
  }
}

template<typename B, typename A>
inplace_update_theta_sketch_alloc<B, A>& inplace_update_theta_sketch_alloc<B, A>::trim() {
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  if (state->num_entries > static_cast<uint32_t>(1 << state->lg_nom_size)) rebuild();
  return *this;
}

template<typename B, typename A>
compact_theta_sketch_alloc<A> inplace_update_theta_sketch_alloc<B, A>::compact(bool ordered, const A& allocator) const {
  // no trimming for now
  std::vector<uint64_t, A> entries(begin(), end(), A(allocator));
  if (ordered) std::sort(entries.begin(), entries.end());
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  return compact_theta_sketch_alloc<A>(state->is_empty, ordered, compute_seed_hash(state->seed), state->theta, std::move(entries));
}

template<typename B, typename A>
auto inplace_update_theta_sketch_alloc<B, A>::begin() -> iterator {
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  return iterator(&state->first_entry, 1 << state->lg_cur_size, 0);
}

template<typename B, typename A>
auto inplace_update_theta_sketch_alloc<B, A>::end() -> iterator {
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  return iterator(nullptr, 0, 1 << state->lg_cur_size);
}

template<typename B, typename A>
auto inplace_update_theta_sketch_alloc<B, A>::begin() const -> const_iterator {
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  return const_iterator(&state->first_entry, 1 << state->lg_cur_size, 0);
}

template<typename B, typename A>
auto inplace_update_theta_sketch_alloc<B, A>::end() const -> const_iterator {
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  return const_iterator(nullptr, 0, 1 << state->lg_cur_size);
}

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::resize() {
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  std::vector<uint64_t, A> entries_copy(begin(), end());
  const uint8_t lg_tgt_size = state->lg_nom_size + 1;
  const uint8_t factor = std::max(1, std::min(static_cast<int>(state->rf), lg_tgt_size - state->lg_cur_size));
  state->lg_cur_size += factor;
  auto state_copy = *state;
  buffer.resize(header_size_bytes() + table_size_bytes(state->lg_cur_size));
  state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  *state = state_copy;
  Entry* entries = &state->first_entry;
  std::fill(entries, entries + (1 << state->lg_cur_size), 0); // might be excessive
  for (auto entry: entries_copy) {
    auto result = Base::find(entries, state->lg_cur_size, entry);
    if (!result.second) *result.first = entry;
  }
}

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::rebuild() {
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(buffer.data());
  const uint32_t nominal_size = 1 << state->lg_nom_size;
  const size_t table_size = 1 << state->lg_cur_size;
  const size_t offset = table_size - state->num_entries;
  Entry* entries = &state->first_entry;
  std::nth_element(entries, entries + nominal_size + offset, entries + table_size);
  state->theta = entries[nominal_size + offset];
  std::vector<uint64_t, A> new_entries(entries + offset, entries + nominal_size + offset);
  std::fill(entries + offset, entries + table_size, 0);
  for (auto entry: new_entries) {
    auto result = Base::find(entries, state->lg_cur_size, entry);
    if (!result.second) *result.first = entry;
  }
  state->num_entries = nominal_size;
}

template<typename B, typename A>
size_t inplace_update_theta_sketch_alloc<B, A>::header_size_bytes() {
  return sizeof(inplace_update_theta_sketch_state) - sizeof(inplace_update_theta_sketch_state::first_entry);
}

template<typename B, typename A>
size_t inplace_update_theta_sketch_alloc<B, A>::table_size_bytes(uint8_t lg_k) {
  return (1 << lg_k) * sizeof(uint64_t);
}

// builder

template<typename B, typename A>
inplace_update_theta_sketch_alloc<B, A>::builder::builder(): theta_base_builder<builder, A>(A()) {}

template<typename B, typename A>
void inplace_update_theta_sketch_alloc<B, A>::builder::initialize(B& buffer) const {
  inplace_update_theta_sketch_alloc::initialize(buffer, this->starting_lg_size(), this->lg_k_, this->rf_, this->starting_theta(), this->seed_);
}

} /* namespace datasketches */

#endif
