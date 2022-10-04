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

#ifndef QUANTILE_SKETCH_SORTED_VIEW_IMPL_HPP_
#define QUANTILE_SKETCH_SORTED_VIEW_IMPL_HPP_

#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace datasketches {

template<typename T, typename C, typename A>
quantile_sketch_sorted_view<T, C, A>::quantile_sketch_sorted_view(uint32_t num, const A& allocator):
total_weight_(0),
entries_(allocator)
{
  entries_.reserve(num);
}

template<typename T, typename C, typename A>
template<typename Iterator>
void quantile_sketch_sorted_view<T, C, A>::add(Iterator first, Iterator last, uint64_t weight) {
  const size_t size_before = entries_.size();
  for (auto it = first; it != last; ++it) entries_.push_back(Entry(ref_helper(*it), weight));
  if (size_before > 0) {
    Container tmp(entries_.get_allocator());
    tmp.reserve(entries_.capacity());
    std::merge(
        entries_.begin(), entries_.begin() + size_before,
        entries_.begin() + size_before, entries_.end(),
        std::back_inserter(tmp), compare_pairs_by_first()
    );
    std::swap(tmp, entries_);
  }
}

template<typename T, typename C, typename A>
void quantile_sketch_sorted_view<T, C, A>::convert_to_cummulative() {
  for (auto& entry: entries_) {
    total_weight_ += entry.second;
    entry.second = total_weight_;
  }
}

template<typename T, typename C, typename A>
double quantile_sketch_sorted_view<T, C, A>::get_rank(const T& item, bool inclusive) const {
  auto it = inclusive ?
      std::upper_bound(entries_.begin(), entries_.end(), Entry(ref_helper(item), 0), compare_pairs_by_first())
    : std::lower_bound(entries_.begin(), entries_.end(), Entry(ref_helper(item), 0), compare_pairs_by_first());
  // we need item just before
  if (it == entries_.begin()) return 0;
  --it;
  return static_cast<double>(it->second) / total_weight_;
}

template<typename T, typename C, typename A>
auto quantile_sketch_sorted_view<T, C, A>::get_quantile(double rank, bool inclusive) const -> quantile_return_type {
  uint64_t weight = inclusive ? std::ceil(rank * total_weight_) : rank * total_weight_;
  auto it = inclusive ?
      std::lower_bound(entries_.begin(), entries_.end(), make_dummy_entry<T>(weight), compare_pairs_by_second())
    : std::upper_bound(entries_.begin(), entries_.end(), make_dummy_entry<T>(weight), compare_pairs_by_second());
  if (it == entries_.end()) return deref_helper(entries_[entries_.size() - 1].first);
  return deref_helper(it->first);
}

template<typename T, typename C, typename A>
auto quantile_sketch_sorted_view<T, C, A>::begin() const -> const_iterator {
  return const_iterator(entries_.begin(), entries_.begin());
}

template<typename T, typename C, typename A>
auto quantile_sketch_sorted_view<T, C, A>::end() const -> const_iterator {
  return const_iterator(entries_.end(), entries_.begin());
}

template<typename T, typename C, typename A>
size_t quantile_sketch_sorted_view<T, C, A>::size() const {
  return entries_.size();
}

} /* namespace datasketches */

#endif
