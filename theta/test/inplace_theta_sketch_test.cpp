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

#include <fstream>

#include <catch.hpp>
#include <inplace_update_theta_sketch.hpp>

#include <theta_union.hpp>

namespace datasketches {

#ifdef TEST_BINARY_INPUT_PATH
const std::string inputPath = TEST_BINARY_INPUT_PATH;
#else
const std::string inputPath = "test/";
#endif

TEST_CASE("inplace theta sketch: empty", "[theta_sketch]") {
  std::vector<char> buf;
  inplace_update_theta_sketch::builder().initialize(buf);
  auto update_sketch = inplace_update_theta_sketch(buf);
  auto compact_sketch = update_sketch.compact();
  REQUIRE(compact_sketch.is_empty());
  REQUIRE_FALSE(compact_sketch.is_estimation_mode());
  REQUIRE(compact_sketch.get_theta() == 1.0);
  REQUIRE(compact_sketch.get_estimate() == 0.0);
  REQUIRE(compact_sketch.get_lower_bound(1) == 0.0);
  REQUIRE(compact_sketch.get_upper_bound(1) == 0.0);
}

TEST_CASE("inplace theta sketch: non empty no retained keys", "[theta_sketch]") {
  std::vector<char> buf;
  inplace_update_theta_sketch::builder().set_p(0.001).initialize(buf);
  auto update_sketch = inplace_update_theta_sketch(buf);
  update_sketch.update(1);
  compact_theta_sketch compact_sketch = update_sketch.compact();
  REQUIRE(compact_sketch.get_num_retained() == 0);
  REQUIRE_FALSE(compact_sketch.is_empty());
  REQUIRE(compact_sketch.is_estimation_mode());
  REQUIRE(compact_sketch.get_estimate() == 0.0);
  REQUIRE(compact_sketch.get_lower_bound(1) == 0.0);
  REQUIRE(compact_sketch.get_upper_bound(1) > 0);
}

TEST_CASE("inplace theta sketch: single item", "[theta_sketch]") {
  std::vector<char> buf;
  inplace_update_theta_sketch::builder().initialize(buf);
  auto update_sketch = inplace_update_theta_sketch(buf);
  update_sketch.update("a");
  compact_theta_sketch compact_sketch = update_sketch.compact();
  REQUIRE_FALSE(compact_sketch.is_empty());
  REQUIRE_FALSE(compact_sketch.is_estimation_mode());
  REQUIRE(compact_sketch.get_theta() == 1.0);
  REQUIRE(compact_sketch.get_estimate() == 1.0);
  REQUIRE(compact_sketch.get_lower_bound(1) == 1.0);
  REQUIRE(compact_sketch.get_upper_bound(1) == 1.0);
}

TEST_CASE("inplace theta sketch: exact", "[theta_sketch]") {
  std::vector<char> buf;
  inplace_update_theta_sketch::builder().initialize(buf);
  auto update_sketch = inplace_update_theta_sketch(buf);
  const int n = 2000;
  for (int i = 0; i < n; i++) update_sketch.update(i);
  auto compact_sketch = update_sketch.compact();
  REQUIRE_FALSE(compact_sketch.is_empty());
  REQUIRE_FALSE(compact_sketch.is_estimation_mode());
  REQUIRE(compact_sketch.get_theta() == 1.0);
  REQUIRE(compact_sketch.get_estimate() == n);
  REQUIRE(compact_sketch.get_lower_bound(1) == n);
  REQUIRE(compact_sketch.get_upper_bound(1) == n);
}

TEST_CASE("inplace theta sketch: estimation", "[theta_sketch]") {
  std::vector<char> buf;
  inplace_update_theta_sketch::builder().initialize(buf);
  auto update_sketch = inplace_update_theta_sketch(buf);
  const int n = 8192;
  for (int i = 0; i < n; ++i) update_sketch.update(i);
  auto compact_sketch = update_sketch.compact();
  REQUIRE_FALSE(compact_sketch.is_empty());
  REQUIRE(compact_sketch.is_ordered());
  REQUIRE(compact_sketch.is_estimation_mode());
  REQUIRE(compact_sketch.get_theta() < 1.0);
  REQUIRE(compact_sketch.get_estimate() == Approx((double) n).margin(n * 0.01));
  REQUIRE(compact_sketch.get_lower_bound(1) < compact_sketch.get_estimate());
  REQUIRE(compact_sketch.get_upper_bound(1) > compact_sketch.get_estimate());
  REQUIRE(compact_sketch.get_num_retained() >= 1 << 12);

  compact_sketch = update_sketch.trim().compact();
  REQUIRE_FALSE(compact_sketch.is_empty());
  REQUIRE(compact_sketch.is_ordered());
  REQUIRE(compact_sketch.is_estimation_mode());
  REQUIRE(compact_sketch.get_theta() < 1.0);
  REQUIRE(compact_sketch.get_estimate() == Approx((double) n).margin(n * 0.01));
  REQUIRE(compact_sketch.get_lower_bound(1) < compact_sketch.get_estimate());
  REQUIRE(compact_sketch.get_upper_bound(1) > compact_sketch.get_estimate());
  REQUIRE(compact_sketch.get_num_retained() == 1 << 12);
}

TEST_CASE("inplace theta sketch: compare with compact sketch in estimation mode from java", "[theta_sketch]") {
  // the same construction process in Java must have produced exactly the same sketch
  std::ifstream is;
  is.exceptions(std::ios::failbit | std::ios::badbit);
  is.open(inputPath + "theta_compact_estimation_from_java.sk", std::ios::binary);
  auto sketch_from_java = compact_theta_sketch::deserialize(is);
  REQUIRE_FALSE(sketch_from_java.is_empty());
  REQUIRE(sketch_from_java.is_estimation_mode());
  REQUIRE(sketch_from_java.is_ordered());
  REQUIRE(sketch_from_java.get_num_retained() == 4342);
  REQUIRE(sketch_from_java.get_theta() == Approx(0.531700444213199).margin(1e-10));
  REQUIRE(sketch_from_java.get_estimate() == Approx(8166.25234614053).margin(1e-10));
  REQUIRE(sketch_from_java.get_lower_bound(2) == Approx(7996.956955317471).margin(1e-10));
  REQUIRE(sketch_from_java.get_upper_bound(2) == Approx(8339.090301078124).margin(1e-10));

  std::vector<char> buf;
  inplace_update_theta_sketch::builder().initialize(buf);
  auto update_sketch = inplace_update_theta_sketch(buf);
  const int n = 8192;
  for (int i = 0; i < n; i++) update_sketch.update(i);
  auto compact_sketch = update_sketch.compact();
  REQUIRE(sketch_from_java.get_num_retained() == compact_sketch.get_num_retained());
  REQUIRE(sketch_from_java.get_theta() == Approx(compact_sketch.get_theta()).margin(1e-10));
  REQUIRE(sketch_from_java.get_estimate() == Approx(compact_sketch.get_estimate()).margin(1e-10));
  REQUIRE(sketch_from_java.get_lower_bound(1) == Approx(compact_sketch.get_lower_bound(1)).margin(1e-10));
  REQUIRE(sketch_from_java.get_upper_bound(1) == Approx(compact_sketch.get_upper_bound(1)).margin(1e-10));
  REQUIRE(sketch_from_java.get_lower_bound(2) == Approx(compact_sketch.get_lower_bound(2)).margin(1e-10));
  REQUIRE(sketch_from_java.get_upper_bound(2) == Approx(compact_sketch.get_upper_bound(2)).margin(1e-10));
  REQUIRE(sketch_from_java.get_lower_bound(3) == Approx(compact_sketch.get_lower_bound(3)).margin(1e-10));
  REQUIRE(sketch_from_java.get_upper_bound(3) == Approx(compact_sketch.get_upper_bound(3)).margin(1e-10));
  // the sketches are ordered, so the iteration sequence must match exactly
  auto iter = sketch_from_java.begin();
  for (const auto& key: compact_sketch) {
    REQUIRE(*iter == key);
    ++iter;
  }
}

TEST_CASE("inplace theta sketch: merge exact disjoint", "[theta_sketch]") {
  std::vector<char> buf1;
  inplace_update_theta_sketch::builder().initialize(buf1);
  auto update_sketch1 = inplace_update_theta_sketch(buf1);
  update_sketch1.update(1);
  update_sketch1.update(2);
  update_sketch1.update(3);
  update_sketch1.update(4);

  std::vector<char> buf2;
  inplace_update_theta_sketch::builder().initialize(buf2);
  auto update_sketch2 = inplace_update_theta_sketch(buf2);
  update_sketch2.update(5);
  update_sketch2.update(6);
  update_sketch2.update(7);
  update_sketch2.update(8);

  update_sketch1.merge(buf2.data());
  auto compact_sketch = update_sketch1.compact();
  REQUIRE_FALSE(compact_sketch.is_empty());
  REQUIRE(compact_sketch.is_ordered());
  REQUIRE_FALSE(compact_sketch.is_estimation_mode());
  REQUIRE(compact_sketch.get_theta() == 1.0);
  REQUIRE(compact_sketch.get_estimate() == 8);
  REQUIRE(compact_sketch.get_lower_bound(1) == 8);
  REQUIRE(compact_sketch.get_upper_bound(1) == 8);
  REQUIRE(compact_sketch.get_num_retained() == 8);
}

TEST_CASE("inplace theta sketch: merge exact half overlap", "[theta_sketch]") {
  std::vector<char> buf1;
  inplace_update_theta_sketch::builder().initialize(buf1);
  auto update_sketch1 = inplace_update_theta_sketch(buf1);
  update_sketch1.update(1);
  update_sketch1.update(2);
  update_sketch1.update(3);
  update_sketch1.update(4);

  std::vector<char> buf2;
  inplace_update_theta_sketch::builder().initialize(buf2);
  auto update_sketch2 = inplace_update_theta_sketch(buf2);
  update_sketch2.update(3);
  update_sketch2.update(4);
  update_sketch2.update(5);
  update_sketch2.update(6);

  update_sketch1.merge(buf2.data());
  auto compact_sketch = update_sketch1.compact();
  REQUIRE_FALSE(compact_sketch.is_empty());
  REQUIRE(compact_sketch.is_ordered());
  REQUIRE_FALSE(compact_sketch.is_estimation_mode());
  REQUIRE(compact_sketch.get_theta() == 1.0);
  REQUIRE(compact_sketch.get_estimate() == 6);
  REQUIRE(compact_sketch.get_lower_bound(1) == 6);
  REQUIRE(compact_sketch.get_upper_bound(1) == 6);
  REQUIRE(compact_sketch.get_num_retained() == 6);
}

TEST_CASE("inplace theta sketch: merge estimation disjoint", "[theta_sketch]") {
  int key = 0;

  std::vector<char> buf1;
  inplace_update_theta_sketch::builder().initialize(buf1);
  auto update_sketch1 = inplace_update_theta_sketch(buf1);
  const int n1 = 8000;
  for (int i = 0; i < n1; ++i) update_sketch1.update(key++);

  std::vector<char> buf2;
  inplace_update_theta_sketch::builder().initialize(buf2);
  auto update_sketch2 = inplace_update_theta_sketch(buf2);
  const int n2 = 16000;
  for (int i = 0; i < n2; ++i) update_sketch2.update(key++);

  update_sketch1.merge(buf2.data());
  const int n = n1 + n2;
  auto result = update_sketch1.compact();

  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_ordered());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_theta() < 1.0);
  REQUIRE(result.get_estimate() == Approx(n).margin(n * 0.02));
  REQUIRE(result.get_lower_bound(1) < result.get_estimate());
  REQUIRE(result.get_upper_bound(1) > result.get_estimate());
  REQUIRE(result.get_num_retained() >= 1 << 12);
}

TEST_CASE("inplace theta sketch: merge compact estimation disjoint", "[theta_sketch]") {
  int key = 0;

  std::vector<char> buf1;
  inplace_update_theta_sketch::builder().initialize(buf1);
  auto update_sketch1 = inplace_update_theta_sketch(buf1);
  const int n1 = 8000;
  for (int i = 0; i < n1; ++i) update_sketch1.update(key++);

  std::vector<char> buf2;
  inplace_update_theta_sketch::builder().initialize(buf2);
  auto update_sketch2 = inplace_update_theta_sketch(buf2);
  const int n2 = 16000;
  for (int i = 0; i < n2; ++i) update_sketch2.update(key++);

  auto bytes = update_sketch2.compact().serialize();
  update_sketch1.merge_compact((char*)bytes.data());
  const int n = n1 + n2;
  auto result = update_sketch1.compact();

  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_ordered());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_theta() < 1.0);
  REQUIRE(result.get_estimate() == Approx(n).margin(n * 0.02));
  REQUIRE(result.get_lower_bound(1) < result.get_estimate());
  REQUIRE(result.get_upper_bound(1) > result.get_estimate());
  REQUIRE(result.get_num_retained() >= 1 << 12);
}

TEST_CASE("inplace theta sketch: merge estimation with overlap", "[theta_sketch]") {
  int key = 0;
  std::vector<char> buf1;
  inplace_update_theta_sketch::builder().initialize(buf1);
  auto update_sketch1 = inplace_update_theta_sketch(buf1);
  const int n1 = 8192;
  for (int i = 0; i < n1; ++i) update_sketch1.update(key++);
  key -= n1 / 2;
  std::vector<char> buf2;
  inplace_update_theta_sketch::builder().initialize(buf2);
  auto update_sketch2 = inplace_update_theta_sketch(buf2);
  const int n2 = 8192;
  for (int i = 0; i < n2; ++i) update_sketch2.update(key++);

  update_sketch1.merge(buf2.data());
  const int n = n1 / 2 + n2;
  auto result = update_sketch1.compact();

  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_ordered());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_theta() < 1.0);
  REQUIRE(result.get_estimate() == Approx(n).margin(n * 0.02));
  REQUIRE(result.get_lower_bound(1) < result.get_estimate());
  REQUIRE(result.get_upper_bound(1) > result.get_estimate());
  REQUIRE(result.get_num_retained() >= 1 << 12);
}

TEST_CASE("inplace theta sketch: merge compact estimation with overlap", "[theta_sketch]") {
  int key = 0;
  std::vector<char> buf1;
  inplace_update_theta_sketch::builder().initialize(buf1);
  auto update_sketch1 = inplace_update_theta_sketch(buf1);
  const int n1 = 8192;
  for (int i = 0; i < n1; ++i) update_sketch1.update(key++);
  key -= n1 / 2;
  std::vector<char> buf2;
  inplace_update_theta_sketch::builder().initialize(buf2);
  auto update_sketch2 = inplace_update_theta_sketch(buf2);
  const int n2 = 8192;
  for (int i = 0; i < n2; ++i) update_sketch2.update(key++);

  auto bytes = update_sketch2.compact().serialize();
  update_sketch1.merge_compact((char*)bytes.data());
  const int n = n1 / 2 + n2;
  auto result = update_sketch1.compact();

  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_ordered());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_theta() < 1.0);
  REQUIRE(result.get_estimate() == Approx(n).margin(n * 0.02));
  REQUIRE(result.get_lower_bound(1) < result.get_estimate());
  REQUIRE(result.get_upper_bound(1) > result.get_estimate());
  REQUIRE(result.get_num_retained() >= 1 << 12);
}

} /* namespace datasketches */
