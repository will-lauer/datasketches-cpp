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

#ifndef REQ_COMMON_HPP_
#define REQ_COMMON_HPP_

#include <random>
#include <chrono>

#include "serde.hpp"
#include "common_defs.hpp"

namespace datasketches {

// TODO: have a common random bit with KLL
static std::independent_bits_engine<std::mt19937, 1, unsigned> req_random_bit(std::chrono::system_clock::now().time_since_epoch().count());

namespace req_constants {
  static const uint16_t MIN_K = 4;
  static const uint8_t INIT_NUM_SECTIONS = 3;
  static const unsigned MULTIPLIER = 2;
}

} /* namespace datasketches */

#endif
