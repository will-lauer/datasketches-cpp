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

namespace datasketches {

template<typename A>
void inplace_update_theta_sketch_alloc<A>::initialize(char* ptr, uint8_t lg_k, float p, uint64_t seed) {
  if (lg_k < theta_constants::MIN_LG_K) {
    throw std::invalid_argument("lg_k must not be less than " + std::to_string(theta_constants::MIN_LG_K) + ": " + std::to_string(lg_k));
  }
  if (lg_k > theta_constants::MAX_LG_K) {
    throw std::invalid_argument("lg_k must not be greater than " + std::to_string(theta_constants::MAX_LG_K) + ": " + std::to_string(lg_k));
  }
  if (p <= 0 || p > 1) throw std::invalid_argument("sampling probability must be between 0 and 1");
  auto state = reinterpret_cast<inplace_update_theta_sketch_state*>(ptr);
  state->is_empty = true;
  state->lg_k = lg_k;
  state->num_entries = 0;
  state->theta = p < 1 ? theta_constants::MAX_THETA * p : theta_constants::MAX_THETA;
  state->seed = seed;
  Entry* entries = &state->first_entry;
  std::fill(entries, entries + (1 << (lg_k + 1)), 0);
}

template<typename A>
inplace_update_theta_sketch_alloc<A>::inplace_update_theta_sketch_alloc(char* ptr):
state(reinterpret_cast<inplace_update_theta_sketch_state*>(ptr)) {}

template<typename A>
void inplace_update_theta_sketch_alloc<A>::update(uint64_t value) {
  update(&value, sizeof(value));
}

template<typename A>
void inplace_update_theta_sketch_alloc<A>::update(const std::string& value) {
  if (value.empty()) return;
  update(value.c_str(), value.length());
}

template<typename A>
size_t inplace_update_theta_sketch_alloc<A>::size_bytes(uint8_t lg_k) {
  return size_u64(lg_k) * sizeof(uint64_t);
}

template<typename A>
void inplace_update_theta_sketch_alloc<A>::update(const void* data, size_t length) {
  state->is_empty = false;
  const uint64_t hash = compute_hash(data, length, state->seed);
  insert_or_ignore(hash);
}

template<typename A>
void inplace_update_theta_sketch_alloc<A>::insert_or_ignore(uint64_t hash) {
  if (hash >= state->theta) return; // hash == 0 is reserved to mark empty slots in the table
  uint64_t* entries = &state->first_entry;
  auto result = Base::find(entries, state->lg_k + 1, hash);
  if (!result.second) {
    *result.first = hash;
    ++state->num_entries;
    if (state->num_entries > Base::get_capacity(state->lg_k + 1, state->lg_k)) {
      rebuild();
    }
  }
}

template<typename A>
void inplace_update_theta_sketch_alloc<A>::merge(const char* ptr) {
  auto other = reinterpret_cast<const inplace_update_theta_sketch_state*>(ptr);
  if (other->is_empty) return;
  // TODO: check seed match
  state->is_empty = false;
  if (state->theta > other->theta) {
    state->theta = other->theta;
    std::vector<uint64_t, A> new_entries;
    for (auto entry: *this) {
      if (entry < state->theta) new_entries.push_back(entry);
    }
    Entry* entries = &state->first_entry;
    std::fill(entries, entries + (1 << (state->lg_k + 1)), 0);
    for (auto entry: new_entries) {
      auto result = Base::find(entries, state->lg_k + 1, entry);
      if (!result.second) *result.first = entry;
    }
    state->num_entries = new_entries.size();
  }
  const Entry* other_entries = &other->first_entry;
  const size_t other_size = 1 << (other->lg_k + 1);
  for (size_t i = 0; i < other_size; ++i) if (other_entries[i]) insert_or_ignore(other_entries[i]);
}

template<typename A>
void inplace_update_theta_sketch_alloc<A>::merge_compact(const char* ptr) {
  // TODO: check serial version and sketch type
  // TODO: check seed match
  if (ptr[COMPACT_SKETCH_FLAGS_BYTE] & (1 << COMPACT_SKETCH_IS_EMPTY_FLAG)) return;
  state->is_empty = false;
  const bool other_has_theta = ptr[COMPACT_SKETCH_PRE_LONGS_BYTE] > 2;
  if (other_has_theta) {
    uint64_t other_theta = reinterpret_cast<const uint64_t*>(ptr)[COMPACT_SKETCH_THETA_U64];
    if (state->theta > other_theta) {
      state->theta = other_theta;
      // TODO: extract common rebuild code
      std::vector<uint64_t, A> new_entries;
      for (auto entry: *this) {
        if (entry < state->theta) new_entries.push_back(entry);
      }
      Entry* entries = &state->first_entry;
      std::fill(entries, entries + (1 << (state->lg_k + 1)), 0);
      for (auto entry: new_entries) {
        auto result = Base::find(entries, state->lg_k + 1, entry);
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

template<typename A>
compact_theta_sketch_alloc<A> inplace_update_theta_sketch_alloc<A>::compact(bool ordered, const A& allocator) const {
  // no trimming for now
  std::vector<uint64_t, A> entries(begin(), end(), A(allocator));
  if (ordered) std::sort(entries.begin(), entries.end());
  return compact_theta_sketch_alloc<A>(state->is_empty, ordered, compute_seed_hash(state->seed), state->theta, std::move(entries));
}

template<typename A>
auto inplace_update_theta_sketch_alloc<A>::begin() -> iterator {
  return iterator(&state->first_entry, 1 << (state->lg_k + 1), 0);
}

template<typename A>
auto inplace_update_theta_sketch_alloc<A>::end() -> iterator {
  return iterator(nullptr, 0, 1 << (state->lg_k + 1));
}

template<typename A>
auto inplace_update_theta_sketch_alloc<A>::begin() const -> const_iterator {
  return const_iterator(&state->first_entry, 1 << (state->lg_k + 1), 0);
}

template<typename A>
auto inplace_update_theta_sketch_alloc<A>::end() const -> const_iterator {
  return const_iterator(nullptr, 0, 1 << (state->lg_k + 1));
}

template<typename A>
void inplace_update_theta_sketch_alloc<A>::rebuild() {
  const uint32_t nominal_size = 1 << state->lg_k;
  const size_t table_size = nominal_size * 2;
  const size_t offset = table_size - state->num_entries;
  Entry* entries = &state->first_entry;
  std::nth_element(entries, entries + nominal_size + offset, entries + table_size);
  state->theta = entries[nominal_size + offset];
  std::vector<uint64_t, A> new_entries(entries + offset, entries + nominal_size + offset);
  std::fill(entries + offset, entries + table_size, 0);
  for (auto entry: new_entries) {
    auto result = Base::find(entries, state->lg_k + 1, entry);
    if (!result.second) *result.first = entry;
  }
  state->num_entries = nominal_size;
}

template<typename A>
size_t inplace_update_theta_sketch_alloc<A>::size_u64(uint8_t lg_k) {
  const size_t header_size_u64 = sizeof(inplace_update_theta_sketch_state) / sizeof(uint64_t) - 1;
  const size_t table_size_u64 = 1 << (lg_k + 1);
  return header_size_u64 + table_size_u64;
}

} /* namespace datasketches */

#endif
