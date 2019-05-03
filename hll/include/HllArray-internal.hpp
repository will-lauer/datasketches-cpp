/*
 * Copyright 2018, Yahoo! Inc. Licensed under the terms of the
 * Apache License 2.0. See LICENSE file at the project root for terms.
 */

#ifndef _HLLARRAY_INTERNAL_HPP_
#define _HLLARRAY_INTERNAL_HPP_

#include "HllArray.hpp"
#include "HllUtil.hpp"
#include "HarmonicNumbers.hpp"
#include "CubicInterpolation.hpp"
#include "CompositeInterpolationXTable.hpp"
//#include "RelativeErrorTables.hpp"
#include "CouponList.hpp"

#include <cstring>
#include <cmath>
#include <stdexcept>
#include <string>

namespace datasketches {

template<typename A>
HllArray<A>::HllArray(const int lgConfigK, const TgtHllType tgtHllType)
  : HllSketchImpl<A>(lgConfigK, tgtHllType, CurMode::HLL) {
  hipAccum = 0.0;
  kxq0 = 1 << lgConfigK;
  kxq1 = 0.0;
  curMin = 0;
  numAtCurMin = 1 << lgConfigK;
  oooFlag = false;
  hllByteArr = nullptr; // allocated in derived class
}

template<typename A>
HllArray<A>::HllArray(const HllArray<A>& that)
  : HllSketchImpl<A>(that.lgConfigK, that.tgtHllType, CurMode::HLL) {
  hipAccum = that.getHipAccum();
  kxq0 = that.getKxQ0();
  kxq1 = that.getKxQ1();
  curMin = that.getCurMin();
  numAtCurMin = that.getNumAtCurMin();
  oooFlag = that.isOutOfOrderFlag();

  // can determine length, so allocate here
  int arrayLen = that.getHllByteArrBytes();
  typedef typename std::allocator_traits<A>::template rebind_alloc<uint8_t> uint8Alloc;
  hllByteArr = uint8Alloc().allocate(arrayLen);
  std::copy(that.hllByteArr, that.hllByteArr + arrayLen, hllByteArr);
}

template<typename A>
HllArray<A>::~HllArray() {
  typedef typename std::allocator_traits<A>::template rebind_alloc<uint8_t> uint8Alloc;
  uint8Alloc().deallocate(hllByteArr, getHllByteArrBytes());
  //delete [] hllByteArr;
}

template<typename A>
HllArray<A>* HllArray<A>::copyAs(const TgtHllType tgtHllType) const {
  if (tgtHllType == this->getTgtHllType()) {
    return static_cast<HllArray*>(copy());
  }
  if (tgtHllType == TgtHllType::HLL_4) {
    //return Conversions::convertToHll4(*this);
    return HllSketchImplFactory<A>::convertToHll4(*this);
  } else if (tgtHllType == TgtHllType::HLL_6) {
    //return Conversions::convertToHll6(*this);
    return HllSketchImplFactory<A>::convertToHll6(*this);
  } else { // tgtHllType == HLL_8
    //return Conversions::convertToHll8(*this);
    return HllSketchImplFactory<A>::convertToHll8(*this);
  }
}

template<typename A>
HllArray<A>* HllArray<A>::newHll(const void* bytes, size_t len) {
  if (len < HllUtil<A>::HLL_BYTE_ARR_START) {
    throw std::invalid_argument("Input data length insufficient to hold HLL array");
  }

  const uint8_t* data = static_cast<const uint8_t*>(bytes);
  if (data[HllUtil<A>::PREAMBLE_INTS_BYTE] != HllUtil<A>::HLL_PREINTS) {
    throw std::invalid_argument("Incorrect number of preInts in input stream");
  }
  if (data[HllUtil<A>::SER_VER_BYTE] != HllUtil<A>::SER_VER) {
    throw std::invalid_argument("Wrong ser ver in input stream");
  }
  if (data[HllUtil<A>::FAMILY_BYTE] != HllUtil<A>::FAMILY_ID) {
    throw std::invalid_argument("Input array is not an HLL sketch");
  }

  CurMode curMode = extractCurMode(data[HllUtil<A>::MODE_BYTE]);
  if (curMode != HLL) {
    throw std::invalid_argument("Calling HLL array construtor with non-HLL mode data");
  }

  TgtHllType tgtHllType = extractTgtHllType(data[HllUtil<A>::MODE_BYTE]);
  bool oooFlag = ((data[HllUtil<A>::FLAGS_BYTE] & HllUtil<A>::OUT_OF_ORDER_FLAG_MASK) ? true : false);
  bool comapctFlag = ((data[HllUtil<A>::FLAGS_BYTE] & HllUtil<A>::COMPACT_FLAG_MASK) ? true : false);

  const int lgK = (int) data[HllUtil<A>::LG_K_BYTE];
  const int curMin = (int) data[HllUtil<A>::HLL_CUR_MIN_BYTE];

  HllArray<A>* sketch = HllSketchImplFactory<A>::newHll(lgK, tgtHllType);
  sketch->putCurMin(curMin);
  sketch->putOutOfOrderFlag(oooFlag);

  int arrayBytes = sketch->getHllByteArrBytes();
  if (len < HllUtil<A>::HLL_BYTE_ARR_START + arrayBytes) {
    throw std::invalid_argument("Input array too small to hold sketch image");
  }

  double hip, kxq0, kxq1;
  std::memcpy(&hip, data + HllUtil<A>::HIP_ACCUM_DOUBLE, sizeof(double));
  std::memcpy(&kxq0, data + HllUtil<A>::KXQ0_DOUBLE, sizeof(double));
  std::memcpy(&kxq1, data + HllUtil<A>::KXQ1_DOUBLE, sizeof(double));
  sketch->putHipAccum(hip);
  sketch->putKxQ0(kxq0);
  sketch->putKxQ1(kxq1);

  int numAtCurMin, auxCount;
  std::memcpy(&numAtCurMin, data + HllUtil<A>::CUR_MIN_COUNT_INT, sizeof(int));
  std::memcpy(&auxCount, data + HllUtil<A>::AUX_COUNT_INT, sizeof(int));
  sketch->putNumAtCurMin(numAtCurMin);

  std::memcpy(sketch->hllByteArr, data + HllUtil<A>::HLL_BYTE_ARR_START, sketch->getHllByteArrBytes());

  if (auxCount > 0) { // necessarily TgtHllType == HLL_4
    int auxLgIntArrSize = (int) data[4];
    const size_t offset = HllUtil<A>::HLL_BYTE_ARR_START + sketch->getHllByteArrBytes();
    const uint8_t* auxDataStart = data + offset;
    AuxHashMap<A>* auxHashMap = AuxHashMap<A>::deserialize(auxDataStart, len - offset, lgK, auxCount, auxLgIntArrSize, comapctFlag);
    ((Hll4Array<A>*)sketch)->putAuxHashMap(auxHashMap);
  }

  return sketch;
}

template<typename A>
HllArray<A>* HllArray<A>::newHll(std::istream& is) {
  uint8_t listHeader[8];
  is.read((char*)listHeader, 8 * sizeof(uint8_t));

  if (listHeader[HllUtil<A>::PREAMBLE_INTS_BYTE] != HllUtil<A>::HLL_PREINTS) {
    throw std::invalid_argument("Incorrect number of preInts in input stream");
  }
  if (listHeader[HllUtil<A>::SER_VER_BYTE] != HllUtil<A>::SER_VER) {
    throw std::invalid_argument("Wrong ser ver in input stream");
  }
  if (listHeader[HllUtil<A>::FAMILY_BYTE] != HllUtil<A>::FAMILY_ID) {
    throw std::invalid_argument("Input stream is not an HLL sketch");
  }

  CurMode curMode = extractCurMode(listHeader[HllUtil<A>::MODE_BYTE]);
  if (curMode != HLL) {
    throw std::invalid_argument("Calling HLL construtor with non-HLL mode data");
  }

  TgtHllType tgtHllType = extractTgtHllType(listHeader[HllUtil<A>::MODE_BYTE]);
  bool oooFlag = ((listHeader[HllUtil<A>::FLAGS_BYTE] & HllUtil<A>::OUT_OF_ORDER_FLAG_MASK) ? true : false);
  bool comapctFlag = ((listHeader[HllUtil<A>::FLAGS_BYTE] & HllUtil<A>::COMPACT_FLAG_MASK) ? true : false);

  const int lgK = (int) listHeader[HllUtil<A>::LG_K_BYTE];
  const int curMin = (int) listHeader[HllUtil<A>::HLL_CUR_MIN_BYTE];

  HllArray* sketch = HllSketchImplFactory<A>::newHll(lgK, tgtHllType);
  sketch->putCurMin(curMin);
  sketch->putOutOfOrderFlag(oooFlag);

  double hip, kxq0, kxq1;
  is.read((char*)&hip, sizeof(hip));
  is.read((char*)&kxq0, sizeof(kxq0));
  is.read((char*)&kxq1, sizeof(kxq1));
  sketch->putHipAccum(hip);
  sketch->putKxQ0(kxq0);
  sketch->putKxQ1(kxq1);

  int numAtCurMin, auxCount;
  is.read((char*)&numAtCurMin, sizeof(numAtCurMin));
  is.read((char*)&auxCount, sizeof(auxCount));
  sketch->putNumAtCurMin(numAtCurMin);
  
  is.read((char*)sketch->hllByteArr, sketch->getHllByteArrBytes());
  
  if (auxCount > 0) { // necessarily TgtHllType == HLL_4
    int auxLgIntArrSize = (int) listHeader[4];
    AuxHashMap<A>* auxHashMap = AuxHashMap<A>::deserialize(is, lgK, auxCount, auxLgIntArrSize, comapctFlag);
    ((Hll4Array<A>*)sketch)->putAuxHashMap(auxHashMap);
  }

  return sketch;
}

template<typename A>
std::pair<std::unique_ptr<uint8_t>, const size_t> HllArray<A>::serialize(bool compact) const {
  const size_t sketchSizeBytes = (compact ? getCompactSerializationBytes() : getUpdatableSerializationBytes());
  typedef typename std::allocator_traits<A>::template rebind_alloc<uint8_t> uint8Alloc;
  std::unique_ptr<uint8_t> byteArr(
    uint8Alloc().allocate(sketchSizeBytes),
    [sketchSizeBytes](uint8_t p){ uint8Alloc().deallocate(p, sketchSizeBytes); }
    );

  uint8_t* bytes = byteArr.get();
  AuxHashMap<A>* auxHashMap = getAuxHashMap();

  bytes[HllUtil<A>::PREAMBLE_INTS_BYTE] = static_cast<uint8_t>(getPreInts());
  bytes[HllUtil<A>::SER_VER_BYTE] = static_cast<uint8_t>(HllUtil<A>::SER_VER);
  bytes[HllUtil<A>::FAMILY_BYTE] = static_cast<uint8_t>(HllUtil<A>::FAMILY_ID);
  bytes[HllUtil<A>::LG_K_BYTE] = static_cast<uint8_t>(this->lgConfigK);
  bytes[HllUtil<A>::LG_ARR_BYTE] = static_cast<uint8_t>(auxHashMap == nullptr ? 0 : auxHashMap->getLgAuxArrInts());
  bytes[HllUtil<A>::FLAGS_BYTE] = this->makeFlagsByte(compact);
  bytes[HllUtil<A>::HLL_CUR_MIN_BYTE] = static_cast<uint8_t>(curMin);
  bytes[HllUtil<A>::MODE_BYTE] = this->makeModeByte();

  std::memcpy(bytes + HllUtil<A>::HIP_ACCUM_DOUBLE, &hipAccum, sizeof(double));
  std::memcpy(bytes + HllUtil<A>::KXQ0_DOUBLE, &kxq0, sizeof(double));
  std::memcpy(bytes + HllUtil<A>::KXQ1_DOUBLE, &kxq1, sizeof(double));
  std::memcpy(bytes + HllUtil<A>::CUR_MIN_COUNT_INT, &numAtCurMin, sizeof(int));
  const int auxCount = (auxHashMap == nullptr ? 0 : auxHashMap->getAuxCount());
  std::memcpy(bytes + HllUtil<A>::AUX_COUNT_INT, &auxCount, sizeof(int));

  const int hllByteArrBytes = getHllByteArrBytes();
  std::memcpy(bytes + getMemDataStart(), hllByteArr, hllByteArrBytes);

  // aux map if HLL_4
  if (this->tgtHllType == HLL_4) {
    bytes += getMemDataStart() + hllByteArrBytes; // start of auxHashMap
    if (auxHashMap != nullptr) {
      if (compact) {
        std::unique_ptr<PairIterator<A>> itr = auxHashMap->getIterator();
        while (itr->nextValid()) {
          const int pairValue = itr->getPair();
          std::memcpy(bytes, &pairValue, sizeof(pairValue));
          bytes += sizeof(pairValue);
        }
      } else {
        std::memcpy(bytes, auxHashMap->getAuxIntArr(), auxHashMap->getUpdatableSizeBytes());
      }
    } else if (!compact) {
      // if updatable, we write even if currently unused so the binary can be wrapped
      int auxBytes = 4 << HllUtil<A>::LG_AUX_ARR_INTS[this->lgConfigK];
      std::fill_n(bytes, auxBytes, 0);
    }
  }

  return std::make_pair(std::move(byteArr), sketchSizeBytes);
}

template<typename A>
void HllArray<A>::serialize(std::ostream& os, const bool compact) const {
  // header
  const uint8_t preInts(getPreInts());
  os.write((char*)&preInts, sizeof(preInts));
  const uint8_t serialVersion(HllUtil<A>::SER_VER);
  os.write((char*)&serialVersion, sizeof(serialVersion));
  const uint8_t familyId(HllUtil<A>::FAMILY_ID);
  os.write((char*)&familyId, sizeof(familyId));
  const uint8_t lgKByte((uint8_t) this->lgConfigK);
  os.write((char*)&lgKByte, sizeof(lgKByte));

  AuxHashMap<A>* auxHashMap = getAuxHashMap();
  uint8_t lgArrByte(0);
  if (auxHashMap != nullptr) {
    lgArrByte = auxHashMap->getLgAuxArrInts();
  }
  os.write((char*)&lgArrByte, sizeof(lgArrByte));

  const uint8_t flagsByte(this->makeFlagsByte(compact));
  os.write((char*)&flagsByte, sizeof(flagsByte));
  const uint8_t curMinByte((uint8_t) curMin);
  os.write((char*)&curMinByte, sizeof(curMinByte));
  const uint8_t modeByte(this->makeModeByte());
  os.write((char*)&modeByte, sizeof(modeByte));

  // estimator data
  os.write((char*)&hipAccum, sizeof(hipAccum));
  os.write((char*)&kxq0, sizeof(kxq0));
  os.write((char*)&kxq1, sizeof(kxq1));

  // array data
  os.write((char*)&numAtCurMin, sizeof(numAtCurMin));

  const int auxCount = (auxHashMap == nullptr ? 0 : auxHashMap->getAuxCount());
  os.write((char*)&auxCount, sizeof(auxCount));
  os.write((char*)hllByteArr, getHllByteArrBytes());

  // aux map if HLL_4
  if (this->tgtHllType == HLL_4) {
    if (auxHashMap != nullptr) {
      if (compact) {
        std::unique_ptr<PairIterator<A>> itr = auxHashMap->getIterator();
        while (itr->nextValid()) {
          const int pairValue = itr->getPair();
          os.write((char*)&pairValue, sizeof(pairValue));
        }
      } else {
        os.write((char*)auxHashMap->getAuxIntArr(), auxHashMap->getUpdatableSizeBytes());
      }
    } else if (!compact) {
      // if updatable, we write even if currently unused so the binary can be wrapped      
      int auxBytes = 4 << HllUtil<A>::LG_AUX_ARR_INTS[this->lgConfigK];
      std::fill_n(std::ostreambuf_iterator<char>(os), auxBytes, 0);
    }
  }
}

template<typename A>
double HllArray<A>::getEstimate() const {
  if (oooFlag) {
    return getCompositeEstimate();
  }
  return getHipAccum();
}

// HLL UPPER AND LOWER BOUNDS

/*
 * The upper and lower bounds are not symmetric and thus are treated slightly differently.
 * For the lower bound, when the unique count is <= k, LB >= numNonZeros, where
 * numNonZeros = k - numAtCurMin AND curMin == 0.
 *
 * For HLL6 and HLL8, curMin is always 0 and numAtCurMin is initialized to k and is decremented
 * down for each valid update until it reaches 0, where it stays. Thus, for these two
 * isomorphs, when numAtCurMin = 0, means the true curMin is > 0 and the unique count must be
 * greater than k.
 *
 * HLL4 always maintains both curMin and numAtCurMin dynamically. Nonetheless, the rules for
 * the very small values <= k where curMin = 0 still apply.
 */
template<typename A>
double HllArray<A>::getLowerBound(const int numStdDev) const {
  HllUtil<A>::checkNumStdDev(numStdDev);
  const int configK = 1 << this->lgConfigK;
  const double numNonZeros = ((curMin == 0) ? (configK - numAtCurMin) : configK);

  double estimate;
  double rseFactor;
  if (oooFlag) {
    estimate = getCompositeEstimate();
    rseFactor = HllUtil<A>::HLL_NON_HIP_RSE_FACTOR;
  } else {
    estimate = hipAccum;
    rseFactor = HllUtil<A>::HLL_HIP_RSE_FACTOR;
  }

  double relErr;
  if (this->lgConfigK > 12) {
    relErr = (numStdDev * rseFactor) / sqrt(configK);
  } else {
    relErr = HllUtil<A>::getRelErr(false, oooFlag, this->lgConfigK, numStdDev);
  }
  return fmax(estimate / (1.0 + relErr), numNonZeros);
}

template<typename A>
double HllArray<A>::getUpperBound(const int numStdDev) const {
  HllUtil<A>::checkNumStdDev(numStdDev);
  const int configK = 1 << this->lgConfigK;

  double estimate;
  double rseFactor;
  if (oooFlag) {
    estimate = getCompositeEstimate();
    rseFactor = HllUtil<A>::HLL_NON_HIP_RSE_FACTOR;
  } else {
    estimate = hipAccum;
    rseFactor = HllUtil<A>::HLL_HIP_RSE_FACTOR;
  }

  double relErr;
  if (this->lgConfigK > 12) {
    relErr = (-1.0) * (numStdDev * rseFactor) / sqrt(configK);
  } else {
    relErr = HllUtil<A>::getRelErr(true, oooFlag, this->lgConfigK, numStdDev);
  }
  return estimate / (1.0 + relErr);
}

/**
 * This is the (non-HIP) estimator.
 * It is called "composite" because multiple estimators are pasted together.
 * @param absHllArr an instance of the AbstractHllArray class.
 * @return the composite estimate
 */
// Original C: again-two-registers.c hhb_get_composite_estimate L1489
template<typename A>
double HllArray<A>::getCompositeEstimate() const {
  const double rawEst = getHllRawEstimate(this->lgConfigK, kxq0 + kxq1);

  const double* xArr = CompositeInterpolationXTable::get_x_arr(this->lgConfigK);
  const int xArrLen = CompositeInterpolationXTable::get_x_arr_length(this->lgConfigK);
  const double yStride = CompositeInterpolationXTable::get_y_stride(this->lgConfigK);

  if (rawEst < xArr[0]) {
    return 0;
  }

  const int xArrLenM1 = xArrLen - 1;

  if (rawEst > xArr[xArrLenM1]) {
    double finalY = yStride * xArrLenM1;
    double factor = finalY / xArr[xArrLenM1];
    return rawEst * factor;
  }

  double adjEst = CubicInterpolation::usingXArrAndYStride(xArr, xArrLen, yStride, rawEst);

  // We need to completely avoid the linear_counting estimator if it might have a crazy value.
  // Empirical evidence suggests that the threshold 3*k will keep us safe if 2^4 <= k <= 2^21.

  if (adjEst > (3 << this->lgConfigK)) { return adjEst; }
  //Alternate call
  //if ((adjEst > (3 << this->lgConfigK)) || ((curMin != 0) || (numAtCurMin == 0)) ) { return adjEst; }

  const double linEst =
      getHllBitMapEstimate(this->lgConfigK, curMin, numAtCurMin);

  // Bias is created when the value of an estimator is compared with a threshold to decide whether
  // to use that estimator or a different one.
  // We conjecture that less bias is created when the average of the two estimators
  // is compared with the threshold. Empirical measurements support this conjecture.

  const double avgEst = (adjEst + linEst) / 2.0;

  // The following constants comes from empirical measurements of the crossover point
  // between the average error of the linear estimator and the adjusted hll estimator
  double crossOver = 0.64;
  if (this->lgConfigK == 4)      { crossOver = 0.718; }
  else if (this->lgConfigK == 5) { crossOver = 0.672; }

  return (avgEst > (crossOver * (1 << this->lgConfigK))) ? adjEst : linEst;
}

template<typename A>
double HllArray<A>::getKxQ0() const {
  return kxq0;
}

template<typename A>
double HllArray<A>::getKxQ1() const {
  return kxq1;
}

template<typename A>
double HllArray<A>::getHipAccum() const {
  return hipAccum;
}

template<typename A>
int HllArray<A>::getCurMin() const {
  return curMin;
}

template<typename A>
int HllArray<A>::getNumAtCurMin() const {
  return numAtCurMin;
}

template<typename A>
void HllArray<A>::putKxQ0(const double kxq0) {
  this->kxq0 = kxq0;
}

template<typename A>
void HllArray<A>::putKxQ1(const double kxq1) {
  this->kxq1 = kxq1;
}

template<typename A>
void HllArray<A>::putHipAccum(const double hipAccum) {
  this->hipAccum = hipAccum;
}

template<typename A>
void HllArray<A>::putCurMin(const int curMin) {
  this->curMin = curMin;
}

template<typename A>
void HllArray<A>::putNumAtCurMin(const int numAtCurMin) {
  this->numAtCurMin = numAtCurMin;
}

template<typename A>
void HllArray<A>::decNumAtCurMin() {
  --numAtCurMin;
}

template<typename A>
void HllArray<A>::addToHipAccum(const double delta) {
  hipAccum += delta;
}

template<typename A>
bool HllArray<A>::isCompact() const {
  return false;
}

template<typename A>
bool HllArray<A>::isEmpty() const {
  const int configK = 1 << this->lgConfigK;
  return (getCurMin() == 0) && (getNumAtCurMin() == configK);
}

template<typename A>
void HllArray<A>::putOutOfOrderFlag(bool flag) {
  oooFlag = flag;
}

template<typename A>
bool HllArray<A>::isOutOfOrderFlag() const {
  return oooFlag;
}

template<typename A>
int HllArray<A>::hll4ArrBytes(const int lgConfigK) {
  return 1 << (lgConfigK - 1);
}

template<typename A>
int HllArray<A>::hll6ArrBytes(const int lgConfigK) {
  const int numSlots = 1 << lgConfigK;
  return ((numSlots * 3) >> 2) + 1;
}

template<typename A>
int HllArray<A>::hll8ArrBytes(const int lgConfigK) {
  return 1 << lgConfigK;
}

template<typename A>
int HllArray<A>::getMemDataStart() const {
  return HllUtil<A>::HLL_BYTE_ARR_START;
}

template<typename A>
int HllArray<A>::getUpdatableSerializationBytes() const {
  return HllUtil<A>::HLL_BYTE_ARR_START + getHllByteArrBytes();
}

template<typename A>
int HllArray<A>::getCompactSerializationBytes() const {
  AuxHashMap<A>* auxHashMap = getAuxHashMap();
  const int auxCountBytes = ((auxHashMap == nullptr) ? 0 : auxHashMap->getCompactSizeBytes());
  return HllUtil<A>::HLL_BYTE_ARR_START + getHllByteArrBytes() + auxCountBytes;
}

template<typename A>
int HllArray<A>::getPreInts() const {
  return HllUtil<A>::HLL_PREINTS;
}

template<typename A>
std::unique_ptr<PairIterator<A>> HllArray<A>::getAuxIterator() const {
  return nullptr;
}

template<typename A>
AuxHashMap<A>* HllArray<A>::getAuxHashMap() const {
  return nullptr;
}

template<typename A>
void HllArray<A>::hipAndKxQIncrementalUpdate(HllArray<A>& host, const int oldValue, const int newValue) {
  if (newValue <= oldValue) {
    throw std::invalid_argument("newValue must be greater than oldValue: " + std::to_string(newValue)
                                + " vs " + std::to_string(oldValue));
  }

  const int configK = 1 << host.getLgConfigK();
  // update hipAccum BEFORE updating kxq0 and kxq1
  double kxq0 = host.getKxQ0();
  double kxq1 = host.getKxQ1();
  host.addToHipAccum(configK / (kxq0 + kxq1));
  // update kxq0 and kxq1; subtract first, then add
  if (oldValue < 32) { host.putKxQ0(kxq0 -= HllUtil<A>::invPow2(oldValue)); }
  else               { host.putKxQ1(kxq1 -= HllUtil<A>::invPow2(oldValue)); }
  if (newValue < 32) { host.putKxQ0(kxq0 += HllUtil<A>::invPow2(newValue)); }
  else               { host.putKxQ1(kxq1 += HllUtil<A>::invPow2(newValue)); }
}

/**
 * Estimator when N is small, roughly less than k log(k).
 * Refer to Wikipedia: Coupon Collector Problem
 * @return the very low range estimate
 */
//In C: again-two-registers.c hhb_get_improved_linear_counting_estimate L1274
template<typename A>
double HllArray<A>::getHllBitMapEstimate(const int lgConfigK, const int curMin, const int numAtCurMin) const {
  const  int configK = 1 << lgConfigK;
  const  int numUnhitBuckets =  ((curMin == 0) ? numAtCurMin : 0);

  //This will eventually go away.
  if (numUnhitBuckets == 0) {
    return configK * log(configK / 0.5);
  }

  const int numHitBuckets = configK - numUnhitBuckets;
  return HarmonicNumbers::getBitMapEstimate(configK, numHitBuckets);
}

//In C: again-two-registers.c hhb_get_raw_estimate L1167
template<typename A>
double HllArray<A>::getHllRawEstimate(const int lgConfigK, const double kxqSum) const {
  const int configK = 1 << lgConfigK;
  double correctionFactor;
  if (lgConfigK == 4) { correctionFactor = 0.673; }
  else if (lgConfigK == 5) { correctionFactor = 0.697; }
  else if (lgConfigK == 6) { correctionFactor = 0.709; }
  else { correctionFactor = 0.7213 / (1.0 + (1.079 / configK)); }
  const double hyperEst = (correctionFactor * configK * configK) / kxqSum;
  return hyperEst;
}

}

#endif // _HLLARRAY_INTERNAL_HPP_
