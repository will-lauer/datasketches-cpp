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
#include <inplace_theta_intersection.hpp>

#include <theta_intersection.hpp>

namespace datasketches {

#ifdef TEST_BINARY_INPUT_PATH
const std::string inputPath = TEST_BINARY_INPUT_PATH;
#else
const std::string inputPath = "test/";
#endif

TEST_CASE("inplace theta intersection: invalid", "[theta_sketch]") {
  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  REQUIRE_THROWS_AS(intersection.get_result(), std::invalid_argument);
}

TEST_CASE("inplace theta intersection compact: empty", "[theta_intersection]") {
  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  update_theta_sketch sketch = update_theta_sketch::builder().build();
  auto bytes = sketch.compact().serialize();
  intersection.intersection_compact(reinterpret_cast<char*>(bytes.data()), bytes.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE(result.get_num_retained() == 0);
  REQUIRE(result.is_empty());
  REQUIRE_FALSE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 0.0);

  intersection.intersection_compact(reinterpret_cast<char*>(bytes.data()), bytes.size());
  result = intersection.get_result();
  REQUIRE(result.get_num_retained() == 0);
  REQUIRE(result.is_empty());
  REQUIRE_FALSE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 0.0);
}

TEST_CASE("inplace theta intersection compact: non empty no retained keys", "[theta_intersection]") {
  update_theta_sketch sketch = update_theta_sketch::builder().set_p(0.001).build();
  sketch.update(1);
  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  auto bytes = sketch.compact().serialize();
  intersection.intersection_compact(reinterpret_cast<char*>(bytes.data()), bytes.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE(result.get_num_retained() == 0);
  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_theta() == Approx(0.001).margin(1e-10));
  REQUIRE(result.get_estimate() == 0.0);

  intersection.intersection_compact(reinterpret_cast<char*>(bytes.data()), bytes.size());
  result = intersection.get_result();
  REQUIRE(result.get_num_retained() == 0);
  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_theta() == Approx(0.001).margin(1e-10));
  REQUIRE(result.get_estimate() == 0.0);
}

TEST_CASE("inplace theta intersection compact: exact mode half overlap unordered", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 1000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact().serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  value = 500;
  for (int i = 0; i < 1000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact().serialize();

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());
  intersection.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());
  compact_theta_sketch result = intersection.get_result();
  //REQUIRE_FALSE(result.is_empty());
  REQUIRE_FALSE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 500.0);
}

TEST_CASE("inplace theta intersection compact: exact mode half overlap ordered", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 1000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact().serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  value = 500;
  for (int i = 0; i < 1000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact().serialize();

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());
  intersection.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE_FALSE(result.is_empty());
  REQUIRE_FALSE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 500.0);
}

TEST_CASE("inplace theta intersection compact: exact mode disjoint unordered", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 1000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact(false).serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  for (int i = 0; i < 1000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact(false).serialize();

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());
  intersection.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE(result.is_empty());
  REQUIRE_FALSE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 0.0);
}

TEST_CASE("inplace theta intersection compact: exact mode disjoint ordered", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 1000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact().serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  for (int i = 0; i < 1000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact().serialize();

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());
  intersection.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE(result.is_empty());
  REQUIRE_FALSE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 0.0);
}

TEST_CASE("inplace theta intersection compact: estimation mode half overlap unordered", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 10000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact(false).serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  value = 5000;
  for (int i = 0; i < 10000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact(false).serialize();

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());
  intersection.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == Approx(5000).margin(5000 * 0.02));
}

TEST_CASE("inplace theta intersection compact: estimation mode half overlap ordered", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 10000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact().serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  value = 5000;
  for (int i = 0; i < 10000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact().serialize();

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());
  intersection.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == Approx(5000).margin(5000 * 0.02));
}

TEST_CASE("inplace theta intersection compact: estimation mode disjoint unordered", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 10000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact(false).serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  for (int i = 0; i < 10000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact(false).serialize();

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());
  intersection.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 0.0);
}

TEST_CASE("inplace theta intersection compact: estimation mode disjoint ordered", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 10000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact().serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  for (int i = 0; i < 10000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact().serialize();

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());
  intersection.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 0.0);
}

TEST_CASE("inplace theta intersection: seed mismatch", "[theta_intersection]") {
  update_theta_sketch sketch = update_theta_sketch::builder().build();
  sketch.update(1); // non-empty should not be ignored
  auto bytes = sketch.compact().serialize();
  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf, 123);
  inplace_theta_intersection intersection(buf);
  REQUIRE_THROWS_AS(intersection.intersection_compact(reinterpret_cast<char*>(bytes.data()), bytes.size()), std::invalid_argument);
}

TEST_CASE("inplace theta intersection: exact mode half overlap", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 1000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact().serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  value = 500;
  for (int i = 0; i < 1000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact().serialize();

  std::vector<char> buf1;
  inplace_theta_intersection::initialize(buf1);
  inplace_theta_intersection intersection1(buf1);
  intersection1.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());

  std::vector<char> buf2;
  inplace_theta_intersection::initialize(buf2);
  inplace_theta_intersection intersection2(buf2);
  intersection2.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection(buf1.data(), buf1.size());
  intersection.intersection(buf2.data(), buf2.size());
  compact_theta_sketch result = intersection.get_result();
  //REQUIRE_FALSE(result.is_empty());
  REQUIRE_FALSE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 500.0);
}

TEST_CASE("inplace theta intersection: exact mode disjoint", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 1000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact(false).serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  for (int i = 0; i < 1000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact(false).serialize();

  std::vector<char> buf1;
  inplace_theta_intersection::initialize(buf1);
  inplace_theta_intersection intersection1(buf1);
  intersection1.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());

  std::vector<char> buf2;
  inplace_theta_intersection::initialize(buf2);
  inplace_theta_intersection intersection2(buf2);
  intersection2.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection(buf1.data(), buf1.size());
  intersection.intersection(buf2.data(), buf2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE(result.is_empty());
  REQUIRE_FALSE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 0.0);
}

TEST_CASE("inplace theta intersection: estimation mode half overlap", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 10000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact(false).serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  value = 5000;
  for (int i = 0; i < 10000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact(false).serialize();

  std::vector<char> buf1;
  inplace_theta_intersection::initialize(buf1);
  inplace_theta_intersection intersection1(buf1);
  intersection1.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());

  std::vector<char> buf2;
  inplace_theta_intersection::initialize(buf2);
  inplace_theta_intersection intersection2(buf2);
  intersection2.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection(buf1.data(), buf1.size());
  intersection.intersection(buf2.data(), buf2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == Approx(5000).margin(5000 * 0.02));
}

TEST_CASE("inplace theta intersection: estimation mode disjoint", "[theta_intersection]") {
  update_theta_sketch sketch1 = update_theta_sketch::builder().build();
  int value = 0;
  for (int i = 0; i < 10000; i++) sketch1.update(value++);
  auto bytes1 = sketch1.compact(false).serialize();

  update_theta_sketch sketch2 = update_theta_sketch::builder().build();
  for (int i = 0; i < 10000; i++) sketch2.update(value++);
  auto bytes2 = sketch2.compact(false).serialize();

  std::vector<char> buf1;
  inplace_theta_intersection::initialize(buf1);
  inplace_theta_intersection intersection1(buf1);
  intersection1.intersection_compact(reinterpret_cast<char*>(bytes1.data()), bytes1.size());

  std::vector<char> buf2;
  inplace_theta_intersection::initialize(buf2);
  inplace_theta_intersection intersection2(buf2);
  intersection2.intersection_compact(reinterpret_cast<char*>(bytes2.data()), bytes2.size());

  std::vector<char> buf;
  inplace_theta_intersection::initialize(buf);
  inplace_theta_intersection intersection(buf);
  intersection.intersection(buf1.data(), buf1.size());
  intersection.intersection(buf2.data(), buf2.size());
  compact_theta_sketch result = intersection.get_result();
  REQUIRE_FALSE(result.is_empty());
  REQUIRE(result.is_estimation_mode());
  REQUIRE(result.get_estimate() == 0.0);
}

} /* namespace datasketches */
